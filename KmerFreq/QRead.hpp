#pragma once

#include <string>
#include <map>
#include <unordered_map>
#include <thread>

class QRead
{
public:
    explicit QRead(std::string filename) : filename_(filename), bad_file_(false)
    {
        cores_ = std::thread::hardware_concurrency(); // indication (CPU cores)
    }
    void kmer_freq(unsigned kmersize, unsigned topcount);
    
    
private:
    unsigned cores_;
    std::string filename_;
    size_t linecount_;
    size_t readlen_;
    bool bad_file_;
    std::once_flag check_file_flag_;
    std::unordered_map<std::string, unsigned> kmers_;
    
    
    void check_file(void);
    void process_read(const std::string &line, unsigned kmersize);
    void generate(unsigned kmersize);
    void generate_paged(unsigned kmersize);
    
};