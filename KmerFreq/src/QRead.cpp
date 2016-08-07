/**
 @file
     QRead.cpp
     Program that produce the most frequent DNA-k-mers of arbitrary length sorted by frequency in a given FASTQ file.
 
 @author
     Name:          Evren KANALICI
     Date:          06/08/2016
     E-Mail:        kanalici.evren@gmail.com
 */
#include <iostream>
#include <algorithm>
#include <sstream>
#include <numeric>
#include <chrono>

#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "QRead.hpp"
#include "Util.hpp"


namespace KmerFreq {

// Flag to enable paged FASTQ file reading
#define PAGING          (1)
// Maxsize size for a sequence-read lenght
#define MAXSIZE_READ    (256)

using namespace std;


/**
    Checks FASTQ file for consistancy;
        - file exists
        - linecount is OK (N x 4)
        - first seq-read is read
 */
void QRead::check_file()
{
    // open file stream
    ifstream infile(filename_);
    if(!infile.good())
    {
        cerr << "! file `" << filename_ << "` does not exists!" << endl;
        bad_file_ = true;
        return;
    }
    
    // get line count
    linecount_ = count(istreambuf_iterator<char>(infile),
                       istreambuf_iterator<char>(), '\n');
    cout << "# line count: " << linecount_ << endl;
    
    // FASTQ file; line count is multiple of 4
    if(linecount_%4 != 0)
    {
        cerr << "! bad FASTQ file!" << endl;
        bad_file_ = true;
        return;
    }
    
    // reset & skip first line (read name)
    infile.seekg(0);
    Util::skip_line(infile);
    
    // Read first seq-read
    string line;
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


/**
    Process & update k-mer occurances for a seq-read line
 
    @param line     - seq-read line
    @param kmersize - kmer (substring) size for occurances
 */
void QRead::process_read(const string &line, unsigned kmersize)
{
    // lock for reentrance
    lock_guard<mutex> kmers_locker(kmers_mutex_);
    
    for(int i=0; i<line.length()-kmersize+1; ++i)
    {
        // get kmer from read
        string kmer = line.substr(i, kmersize);
        // get occurance
        auto slot = kmers_.find(kmer);
        
        // update occurance
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


/**
    Calcuates range offset for n-partition over file map of file_len size
    calculation happen in following fashion;
 
        for n partition, there will be (n+1) offsets [start:end) each starting in FASTQ file block
 
    @param map          - file map
    @param file_len     - total file lenght
    @param n            - partition count
 
    @return             - vector of (n+1) offsets
 
 */
vector<size_t> QRead::qoffset_ranges(char *map, size_t file_len, unsigned n) const
{
    vector<size_t> offsets;
    // first offset
    offsets.push_back(0);
    
    // calcuate in-between ofsets for n
    for(int i=0; i<n-1; ++i)
    {
        size_t offset = (i+1)*(file_len/n);
        // seek till beginnin of the block
        while(map[offset] != '@')
            offset -= 1;
        // in-between offset
        offsets.push_back(offset);
    }
    //last offset
    offsets.push_back(file_len);
    
    return std::move(offsets);
}

/**
    Task for processing ranges of FASTQ file paged-map
    ranges [start:end) calcuated according to n partition #see `qoffset_ranges`
    
    @param map      - file map
    @param start    - start offset of the range
    @param end      - end offset of the range
    @param kmersize - kmer size for processing lines
 
 */
int QRead::qrange_task(char *map, size_t start, size_t end, unsigned kmersize)
{
    unsigned long offset = 0;  // map offset
    char block_offset = 0;  // FASTQ file block(4 lines) offset
    int line_offset = 0;  // in-line offset
    unsigned long cur = 0;  // current char
    char buf_line[MAXSIZE_READ+1] = {'\0'};  // line buffer
    
    // iterate map over the range
    for(offset = start; offset < end; ++offset)
    {
        cur = map[offset];
        
        // end of the line, update offsets
        if(cur == '\n')
        {
            block_offset = (block_offset + 1) % 4;
            
            // finalize line buffer
            buf_line[max(MAXSIZE_READ, line_offset)] = '\0';
            line_offset = -1;
        }
        else
        {
            line_offset++;
        }
        
        // sequcne-read line
        if(block_offset == 1)
        {
            /*
            if(line_offset == 0)
            {
                memset(buf_line, '\0', (MAXSIZE_READ+1));
            }
            */
             
            
            //TODO: switch to batch'd realloc()!
            if(line_offset >= MAXSIZE_READ)
            {
                stringstream ss;
                ss << "[" << this_thread::get_id() << "] " << "MAXSIZE_READ reached!";
                
                throw out_of_range(ss.str());
            }
            
            // update line buffer
            buf_line[line_offset] = (char)cur;
            
        }
        // sequence read
        else if(block_offset == 2 && line_offset == 0)
        {
            // process read line
            string line(buf_line);
            process_read(line, kmersize);
        }
        else {
            //TODO: determistic offset skip? (no guarantee for read lengths)
        }
    }
    
    return EXIT_SUCCESS;
}

/**
    Generates kmer occurances of kmerize using `process_read` for entire FASTQ file
 
    @param kmersize     - kmer size for processing lines
 */
void QRead::generate(unsigned kmersize)
{
    string line;  // seq-read line
    long idx_read = 0;  // read index
    ifstream infile(filename_);  // FASTQ file stream
    
    // reset & skip first line (read name)
    infile.seekg(0);
    Util::skip_line(infile);
    
    // iterate lines in block basis
    while(getline(infile, line))
    {
#if DEBUG
        cout << "+ read-" << ++idx_read << ": ";
        cout << line << endl;
#endif
        // process read line
        process_read(line, kmersize);
        
        // skip block
        Util::skip_line(infile, 3);
    }
}

/**
    Generates kmer occurances of kmerize using `process_read` for entire FASTQ file
    in a concurrent fashion using mem. paging
 
    @param kmersize     - kmer size for processing lines
 */
void QRead::generate_paged(unsigned kmersize)
{
    struct stat sb; // file stat
    int fd;  // file handle
    
    // open file & report error
    if ((fd = open(filename_.c_str(), O_RDONLY)) == -1)
    {
        fatal_error("error opening file!");
    }
    if (fstat(fd, &sb) == -1)
    {
        close(fd);
        fatal_error("error retrieving file!");
    }
    
    // get file lenght
    unsigned long file_len = (unsigned long)lseek(fd, 0, SEEK_END);
    
    // get page map for file & report error
    char *map = (char*)mmap(NULL, file_len, PROT_READ, MAP_SHARED, fd, 0);
    if (map == MAP_FAILED) {
        close(fd);
        fatal_error("error mmapping the file!");
    }
    
    // calcuate offset ranges for tasts
    auto ranges = qoffset_ranges(map, file_len, cores_);
    range_threads_.clear();
    
    // generate range tasts
    for(auto it = ranges.begin(); it != ranges.end()-1; ++it)
    {
        // offset range
        size_t start = *it;
        size_t end = *(it+1);
        
        /*
         * range tasks with ranges to read FASTQ file concurrently
         */
        future<int> fu = async(launch::async, &QRead::qrange_task, this, map, start, end, kmersize);
        range_threads_.emplace_back(std::move(fu));
    }
    
    try
    {
        // wait for range tasts
        for (auto& fu : range_threads_) {
            fu.get();  // ### std::future trick to propagate exception ###
        }
    }
    catch(const exception& e)
    {
        munmap(map, file_len);
        close(fd);
        
        fatal_error(e.what());
    }
    
    
    // unmap & close handle
    munmap(map, file_len);
    close(fd);
    
}


#pragma mark - Public

/**
    Algorith entrance function for calculation & reporting `topcount` frequent k-mers of `kmersize` for FASTQ file format
 
    @param kmerize      - kmer size to process lines
    @param topcount     - top frequent kmers to calculate&report
 */
void QRead::kmer_freq(unsigned kmersize, unsigned topcount)
{
    // check file for consistency & format
    cout << "# checking file: `" << filename_ << "`..." << endl;
    call_once(check_file_flag_, &QRead::check_file, this);
    if(bad_file_) return;
    
    // check kmersize passed
    if(kmersize > readlen_)
    {
        cerr << "! k-mer size is too big!";
        return;
    }
    
    // clear kmer_ occurances data
    kmers_.clear();
    
    // start algo. timer
    cout << "# running algorithm..." << endl;
    auto wcts = std::chrono::system_clock::now();
    
    /*
     * Generate k-mer occurances data
     */
#if PAGING
    generate_paged(kmersize);
#else
    generate(kmersize);
#endif

    // flip occurances map to sort for values
    auto bi_kmers = Util::flip_map(kmers_);
    
    // report algo. timer
    chrono::duration<double> wctduration = (chrono::system_clock::now() - wcts);
    cout << "# time spent: " << wctduration.count() << " secs." << endl;
    
    /*
     * Calcuate total occurances
     * read length can be different (not commonly), so
     * linecount_/4 * (readlenght-kmersize+1) not ALWAYS correct
     */
    const size_t occur_sum = accumulate(kmers_.begin(), kmers_.end(), 0, [](const size_t acc, const std::pair<std::string,size_t>& p) {
        return acc + p.second;
    });
    cout << "# total kmer count: " << occur_sum << endl;
    
    
    // report topcount frequent k-mers
    cout << endl << "# showing top " << topcount << " frequent " << kmersize << "-mers;" << endl;
    cout << "------------------------------------------------------" << endl;
    int idx = 0;
    cout.precision(2);
    cout.setf(ios::fixed);
    for(const auto& kv: bi_kmers)
    {
        cout << kv.second << " :\t %" << 100*kv.first/static_cast<double>(occur_sum) << " - " << kv.first << endl;
        if(++idx == topcount)
        {
            break;
        }
    }
}
 
/* end namespace KmerFreq */
};
