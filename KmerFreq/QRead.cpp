#include <iostream>
#include <algorithm>
#include <numeric>
#include <chrono>

#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "QRead.hpp"
#include "Util.hpp"


using namespace std;

#define PAGING          (1)
#define MAXSIZE_READ    (256)



void QRead::check_file()
{
    ifstream infile(filename_);
    if(!infile.good())
    {
        cerr << "! file `" << filename_ << "` does not exists!" << endl;
        bad_file_ = true;
        return;
    }
    
    
    linecount_ = count(istreambuf_iterator<char>(infile),
                       istreambuf_iterator<char>(), '\n');
    
    
    cout << "# line count: " << linecount_ << endl;
    
    if(linecount_%4 != 0)
    {
        cerr << "! bad FASTQ file!" << endl;
        bad_file_ = true;
        return;
    }
    
    //read first read seq.
    string line;
    infile.seekg(0);
    Util::skip_line(infile);
    
    
    if(getline(infile, line))
    {
        readlen_ = line.length();
        cout << "# read lenght: " << readlen_ << endl;
    }
    else
    {
        cerr << "! bad FASTQ format!" << endl;
        bad_file_ = true;
        return;
    }
    
    cout << "# [file OK]" << endl;
}


void QRead::process_read(const string &line, unsigned kmersize)
{
    lock_guard<mutex> kmers_locker(kmers_mutex_);
    
    for(int i=0; i<line.length()-kmersize+1; ++i)
    {
        string kmer = line.substr(i, kmersize);
        auto slot = kmers_.find(kmer);
        if(slot == kmers_.end())
        {
            kmers_[kmer] = 1;
        }
        else
        {
            slot->second += 1;
        }
    }
}


//(n+1) offsets
vector<size_t> QRead::qoffset_ranges(char *map, size_t file_len, unsigned n)
{
    vector<size_t> offsets;
    offsets.push_back(0);
    for(int i=0; i<n-1; ++i)
    {
        size_t offset = (i+1)*(file_len/n);
        while(map[offset] != '@')
            offset -= 1;
        offsets.push_back(offset);
    }
    offsets.push_back(file_len);
    
    return std::move(offsets);
}

void QRead::qrange_task(char *map, size_t start, size_t end, unsigned kmersize)
{
    unsigned long offset = 0;
    unsigned long cur = 0;
    char buf_line[MAXSIZE_READ+1];
    
    
    char block_offset = 0;
    int line_offset = 0;
    for(offset = start; offset < end; ++offset)
    {
        cur = map[offset];
        if(cur == '\n')
        {
            block_offset = (block_offset + 1) % 4;
            
            buf_line[max(MAXSIZE_READ, line_offset)] = '\0';
            line_offset = -1;
        }
        else
        {
            line_offset++;
        }
        
        if(block_offset == 1)
        {
            
            if(line_offset == 0)
            {
                //memset(buf_line, '\0', (MAXSIZE_READ+1));
            }
            
            if(line_offset > MAXSIZE_READ)
            {
                handle_error("MAXSIZE_READ reached!");
            }
            buf_line[line_offset] = (char)cur;
            
        }
        else if(block_offset == 2 && line_offset == 0) // sequence read
        {
            
            string line(buf_line);
            process_read(line, kmersize);
        }
        else {
            //TODO: determistic skip? (no guarantee for read lenghts)
        }
    }
    
}


void QRead::generate(unsigned kmersize)
{
    string line;
    long idx_read = 0;
    ifstream infile(filename_);
    
    //reset & skip first line
    infile.seekg(0);
    Util::skip_line(infile);
    
    while(getline(infile, line))
    {
#if DEBUG
        cout << "+ read-" << ++idx_read << ": ";
        cout << line << endl;
#endif
        process_read(line, kmersize);
        
        Util::skip_line(infile, 3);
    }
}


void QRead::generate_paged(unsigned kmersize)
{
    
    struct stat sb;
    int fd = open(filename_.c_str(), O_RDONLY);
    
    if (fd == -1)
    {
        handle_error("error opening file!");
    }
    if (fstat(fd, &sb) == -1)
    {
        close(fd);
        handle_error("error retrieving file!");
    }
    
    unsigned long file_len = (unsigned long)lseek(fd, 0, SEEK_END);
    
    //page map
    char *map = (char*)mmap(NULL, file_len, PROT_READ, MAP_SHARED, fd, 0);
    if (map == MAP_FAILED) {
        close(fd);
        handle_error("Error mmapping the file");
    }
    
    
    auto ranges = qoffset_ranges(map, file_len, cores_);
    vector<thread> tasks;
    for(auto it = ranges.begin(); it != ranges.end()-1; ++it)
    {
        size_t start = *it;
        size_t end = *(it+1);
        
        tasks.emplace_back(&QRead::qrange_task, this, map, start, end, kmersize);
    }
    for (auto& t : tasks) {
        t.join();
    }
    
    munmap(map, file_len);
    close(fd);
    

}


#pragma mark - Public


void QRead::kmer_freq(unsigned kmersize, unsigned topcount)
{
    cout << "# checking file: `" << filename_ << "`..." << endl;
    
    call_once(check_file_flag_, &QRead::check_file, this);
    if(bad_file_) return;
    
    
    if(kmersize > readlen_)
    {
        cerr << "! k-mer size is too big!";
        return;
    }
    
    kmers_.clear();
    
    cout << "# running..." << endl;
    auto wcts = std::chrono::system_clock::now();
    
#if PAGING
    generate_paged(kmersize);
#else
    generate(kmersize);
#endif

    //flip occurances
    auto dst = Util::flip_map(kmers_);
    
    chrono::duration<double> wctduration = (chrono::system_clock::now() - wcts);
    cout << "# time spent: " << wctduration.count() << " secs." << endl;
    
    //read length can be different (not commonly), so
    //linecount_/4 * (readlenght-kmersize+1)
    const size_t occur_sum = accumulate(kmers_.begin(), kmers_.end(), 0, [](const size_t acc, const std::pair<std::string,size_t>& p) {
        return acc + p.second;
    });
    cout << "# total kmer count: " << occur_sum << endl;
    cout << endl << "# showing top " << topcount << " frequent kmer;" << endl;
    cout << "-------------------------------------------------------------" << endl;
    
    int idx = 0;
    cout.precision(2);
    cout.setf(ios::fixed);
    for(const auto& kv: dst)
    {
        cout << kv.second << " :\t %" << 100*kv.first/static_cast<double>(occur_sum) << " - " << kv.first << endl;
        if(++idx == topcount)
        {
            break;
        }
    }
    
}
