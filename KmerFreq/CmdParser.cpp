/**
 @file
     CmdParser.cpp
     Program that produce the most frequent DNA-k-mers of arbitrary length sorted by frequency in a given FASTQ file.
 
 @author
     Name:          Evren KANALICI
     Date:          06/08/2016
     E-Mail:        kanalici.evren@gmail.com
 */
#include <sstream>
#include "CmdParser.hpp"

using namespace std;

namespace KmerFreq {


const std::string CmdParser::ARG_FILENAME = "filename";
const std::string CmdParser::ARG_KMERSIZE = "kmersize";
const std::string CmdParser::ARG_TOPCOUNT = "topcount";
const std::string CmdParser::USAGE = "usage : kmer_freq --filename <path> --kmersize <size> --topcount <count>";


/**
    Allocate CmdParser with main() arguments
    
    @param argc - number of args
    @param argv - arg buffer
    
    @throws CmdParserException on invalid argument usage
 
 */
CmdParser::CmdParser(int argc, char** argv) throw(CmdParserException&)
{
    exec_name_ = string(argv[0]);
    
    // valid usage without label
    if(argc == 3+1)
    {
        args_[ARG_FILENAME] = string(argv[1]);
        args_[ARG_KMERSIZE] = string(argv[2]);
        args_[ARG_TOPCOUNT] = string(argv[3]);
    }
    // valid usage with labels
    else if(argc == 6+1)
    {
        string a;
        for(int i=0; i<3; ++i)
        {
            string key = string(argv[i*2+1]);
            _trim_prefix(key);  // trim dashes from command labels
            args_[key] = string(argv[i*2+2]);
        }
    }
    else
    {
        throw CmdParserException(CmdParser::USAGE.c_str());
    }
    
    // validate command arguments
    _validate();
}

/**
    Trims dashes from command argument labels
    
    @param label    - argument label to trim
 
    @throws CmdParserException on bad label
 */
void CmdParser::_trim_prefix(string& label) const throw(CmdParserException&)
{
    if(label.length() == 0) return;
    if(label.at(0) != '-') return;
    
    size_t begin = label.find_first_not_of("-");
    
    if(begin != 2)
    {
        stringstream ss;
        ss << "invalid label: `" << "`" << endl;
        ss << CmdParser::USAGE << endl;
        
        throw CmdParserException(ss.str().c_str());
    }
    
    label = label.substr(begin);
}

/**
    Validates parsed command arguments
    Args should have a min. subset of pre-defined arguments
 */
void CmdParser::_validate() const
{
    if(args_.find(ARG_FILENAME) == args_.end() ||
       args_.find(ARG_KMERSIZE) == args_.end() ||
       args_.find(ARG_TOPCOUNT) == args_.end())
    {

        throw CmdParserException(CmdParser::USAGE.c_str());
    }
}

/**
    Gets parsed arguments map
 
    @return     - arguments map
*/
const std::map<std::string, std::string>& CmdParser::get() const
{
    return args_;
}
    
/* end namespace KmerFreq */
};
