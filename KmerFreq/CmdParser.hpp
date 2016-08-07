#pragma once

#include <string>
#include <map>
#include <exception>

class CmdParserException : public std::exception
{
public:
    CmdParserException(const char* err_message) : err_message_(err_message){}
    
    const char* what() const throw() { return err_message_; }
    
private:
    const char* err_message_;
};


class CmdParser
{
public:
    static const std::string ARG_FILENAME;
    static const std::string ARG_KMERSIZE;
    static const std::string ARG_TOPCOUNT;
    static const std::string USAGE;
    
    explicit CmdParser(int argc, char** argv) throw(CmdParserException&);
    const std::map<std::string, std::string>& get() const;

private:
    std::map<std::string, std::string> args_;
    std::string exec_name_;
    
    void _trim_prefix(std::string& label) const throw(CmdParserException&);
    void _validate() const;
};


