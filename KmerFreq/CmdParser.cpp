#include <sstream>
#include "CmdParser.hpp"

using namespace std;

const std::string CmdParser::ARG_FILENAME = "filename";
const std::string CmdParser::ARG_KMERSIZE = "kmersize";
const std::string CmdParser::ARG_TOPCOUNT = "topcount";
const std::string CmdParser::USAGE = "usage : kmer_freq --filename <path> --kmersize <size> --topcount <count>"; //TODO

CmdParser::CmdParser(int argc, char** argv) throw(CmdParserException&)
{
    exec_name_ = string(argv[0]);
    
    if(argc == 3+1)
    {
        args_[ARG_FILENAME] = string(argv[1]);
        args_[ARG_KMERSIZE] = string(argv[2]);
        args_[ARG_TOPCOUNT] = string(argv[3]);
    }
    else if(argc == 6+1)
    {
        string a;
        for(int i=0; i<3; ++i)
        {
            string key = string(argv[i*2+1]);
            _trim_prefix(key);
            args_[key] = string(argv[i*2+2]);
        }
    }
    else
    {
        throw CmdParserException(CmdParser::USAGE.c_str());
    }
    
    _validate();
}


void CmdParser::_trim_prefix(string& label) const throw(CmdParserException&)
{
    if(label.length() == 0) return;
    if(label.at(0) != '-') return;
    
    size_t begin = label.find_first_not_of("-");
    if (begin == string::npos)
        return;
    
    if(begin != 2)
    {
        stringstream ss;
        ss << "invalid label: `" << "`" << endl;
        ss << CmdParser::USAGE << endl;
        
        throw CmdParserException(ss.str().c_str());
    }
    
    label = label.substr(begin);
}

void CmdParser::_validate() const
{
    if(args_.find(ARG_FILENAME) == args_.end() ||
       args_.find(ARG_KMERSIZE) == args_.end() ||
       args_.find(ARG_TOPCOUNT) == args_.end())
    {
        stringstream ss;
        throw CmdParserException(CmdParser::USAGE.c_str());
    }
}

const std::map<std::string, std::string>& CmdParser::get() const
{
    return args_;
}
