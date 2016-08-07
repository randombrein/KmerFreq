#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <chrono>

#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "QRead.hpp"
#include "Util.hpp"


using namespace std;

#define MAXSIZE_READ    (256)
#define PAGING          (1)
#define handle_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while (0)



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
        handle_error("error opening file!");
    if (fstat(fd, &sb) == -1)           
        handle_error("error retrieving file!");
    
    unsigned long file_len = (unsigned long)lseek(fd, 0, SEEK_END);
    
    
    //page map
    char *map = (char*)mmap(NULL, file_len, PROT_READ, MAP_SHARED, fd, 0);
    if (map == MAP_FAILED) {
        close(fd);
        handle_error("Error mmapping the file");
    }
    
    unsigned long offset = 0;
    unsigned long cur = 0;
    char buf_line[MAXSIZE_READ+1];
    
    // We only care about the 2nd line of each 4-line block, so some bookkeeping is in order.
    char block_offset = 0;
    short line_offset = 0;
    int idx_read = 0;
    for(offset = 0; offset < file_len; ++offset)
    {
        cur = map[offset];
        if(cur == '\n')
        {
            block_offset = (block_offset + 1) % 4;
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
                memset(buf_line, '\0', (MAXSIZE_READ+1));
            }
            
            if(line_offset > MAXSIZE_READ)
            {
                close(fd);
                munmap(map, file_len);
                handle_error("MAXSIZE_READ reached!");
            }
            buf_line[line_offset] = (char)cur;
            
        }
        else if(block_offset == 2 && line_offset == 0) // sequence read
        {
            string line(buf_line);
#if DEBUG
            cout << "+ read-" << ++idx_read << ": ";
            cout << line << endl;
#endif
            process_read(line, kmersize);
        }
        else {
            //TODO: determistic skip?
        }
    }
    close(fd);
    munmap(map, file_len);

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
