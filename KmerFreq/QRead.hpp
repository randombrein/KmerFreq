/**
 @file
     QRead.hpp
     Program that produce the most frequent DNA-k-mers of arbitrary length sorted by frequency in a given FASTQ file.
 
 @author
     Name:          Evren KANALICI
     Date:          06/08/2016
     E-Mail:        kanalici.evren@gmail.com
 */
#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <thread>

namespace KmerFreq {

class QRead
{
public:
    explicit QRead(std::string filename) : filename_(filename), bad_file_(false)
    {
        cores_ = std::thread::hardware_concurrency(); // indication (CPU cores)
    }
    void kmer_freq(unsigned kmersize, unsigned topcount);
    
    
private:
    unsigned cores_;  // # of CPU cores available on the platform
    std::string filename_;  // filename of the FASTQ file
    size_t linecount_;  // linecount of the FASTQ file
    size_t readlen_;  // sequence-read lenght of the FASTQ file
    bool bad_file_;  // file is not good evaluate
    std::once_flag check_file_flag_;  // check file once_flag
    std::mutex kmers_mutex_;  // kmers_ map mutex
    std::unordered_map<std::string, unsigned> kmers_;  // k-mers occurances map
    
    
    void check_file(void);
    void process_read(const std::string &line, unsigned kmersize);
    std::vector<size_t> qoffset_ranges(char *map, size_t file_len, unsigned n) const;
    void qrange_task(char *map, size_t start, size_t end, unsigned kmersize);
    void generate(unsigned kmersize);
    void generate_paged(unsigned kmersize);
    
};

/* end namespace KmerFreq */
};