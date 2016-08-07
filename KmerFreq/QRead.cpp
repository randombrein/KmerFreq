#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <chrono>
#include "QRead.hpp"
#include "Util.hpp"

using namespace std;


size_t count_line(istream &is)
{
    // skip when bad
    if( is.bad() ) return 0;
    // save state
    std::istream::iostate state_backup = is.rdstate();
    // clear state
    is.clear();
    std::istream::streampos pos_backup = is.tellg();
    
    is.seekg(0);
    size_t line_cnt;
    size_t lf_cnt = std::count(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>(), '\n');
    line_cnt = lf_cnt;
    // if the file is not end with '\n' , then line_cnt should plus 1
    is.unget();
    if( is.get() != '\n' ) { ++line_cnt ; }
    
    // recover state
    is.clear() ; // previous reading may set eofbit
    is.seekg(pos_backup);
    is.setstate(state_backup);
    
    return line_cnt;
}

void QRead::skip_line(ifstream& fs, int num=1)
{
    for(int i=0; i<num; ++i)
        fs.ignore(numeric_limits<streamsize>::max(), '\n' );
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

void QRead::check_file()
{
    ifstream infile(filename_);
    if(!infile.good())
    {
        cerr << "file `" << filename_ << "` does not exists!" << endl;
        bad_file_ = true;
        return;
    }
    
    
    linecount_ = count(istreambuf_iterator<char>(infile),
                       istreambuf_iterator<char>(), '\n');
    
    
    cout << "line count: " << linecount_ << endl;
    
    if(linecount_%4 != 0)
    {
        cerr << "bad FASTQ file!" << endl;
        bad_file_ = true;
        return;
    }
    
    //read first read seq.
    string line;
    infile.seekg(0);
    skip_line(infile);
    
    
    if(getline(infile, line))
    {
        readlen_ = line.length();
        cout << "read lenght: " << readlen_ << endl;
    }
    else
    {
        cerr << "bad FASTQ format!" << endl;
        bad_file_ = true;
        return;
    }
    
    cout << "[file OK]" << endl << endl;;
}



void QRead::kmer_freq(unsigned kmersize, unsigned topcount)
{
    string line;
    long idx = 0;
    ifstream infile(filename_);
    once_flag oflag;
    
    if(bad_file_) return;
    
    cout << "reading file: `" << filename_ << "`..." << endl;
    call_once(oflag, &QRead::check_file, this);
    
    
    if(kmersize > readlen_)
    {
        cerr << "kmer size is too big!";
        return;
    }
    
    //reset & skip first line
    infile.seekg(0);
    skip_line(infile);
    
    kmers_.clear();
    
    //...
    
    auto wcts = std::chrono::system_clock::now();
    
    while(getline(infile, line))
    {
        cout << "read-" << ++idx << ": ";
        cout << line << endl << endl;
        process_read(line, kmersize);
        
        skip_line(infile, 3);
    }
    
    
    chrono::duration<double> wctduration = (chrono::system_clock::now() - wcts);
    cout << "time spent: " << wctduration.count() << " seconds [Wall Clock]" << endl;
    
    //...
    
    
    //flip occurances
    cout << endl << endl;
    auto dst = flip_map(kmers_);
    
    //read length can be different (not commonly), so
    //linecount_/4 * (readlenght-kmersize+1)
    const size_t occur_sum = accumulate(kmers_.begin(), kmers_.end(), 0, [](const size_t acc, const std::pair<std::string,size_t>& p) {
        return acc + p.second;
    });
    cout << "occur_sum: " << occur_sum << endl;
    cout << "showing top " << topcount << " frequent occurances" << endl << endl;
    
    idx = 0;
    cout.precision(2);
    cout.setf(ios::fixed);
    for(const auto& kv: dst)
    {
        cout << kv.second << " :\t\t %" << 100*kv.first/static_cast<double>(occur_sum) << " - " << kv.first << endl;
        if(++idx == topcount)
        {
            break;
        }
    }
    
}
