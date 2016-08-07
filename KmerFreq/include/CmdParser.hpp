/**
 @file
     CmdParser.hpp
     Program that produce the most frequent DNA-k-mers of arbitrary length sorted by frequency in a given FASTQ file.
 
 @author
     Name:          Evren KANALICI
     Date:          06/08/2016
     E-Mail:        kanalici.evren@gmail.com
 */
#pragma once

#include <string>
#include <map>
#include <exception>

namespace KmerFreq {

    
class CmdParserException : public std::exception
{
public:
    CmdParserException(const std::string& err_message) : err_message_(err_message.c_str()){}
    CmdParserException(const char* err_message) : err_message_(err_message){}
    
    const char* what() const throw() { return err_message_; }
    
private:
    const char* err_message_;
};


class CmdParser
{
public:
    // Command parser arguments
    static const std::string ARG_FILENAME;
    static const std::string ARG_KMERSIZE;
    static const std::string ARG_TOPCOUNT;
    // Short command usage string
    static const std::string USAGE;
    
    explicit CmdParser(int argc, char** argv) throw(CmdParserException&);
    const std::map<std::string, std::string>& get() const;

private:
    std::map<std::string, std::string> args_;  // command arguments
    std::string exec_name_;  //executable name
    
    void _trim_prefix(std::string& label) const throw(CmdParserException&);
    void _validate() const;
};
  
    
/* end namespace KmerFreq */
};


