#include <iostream>
#include <string>
#include <exception>

#include "CmdParser.hpp"
#include "QRead.hpp"

using namespace std;

int main(int argc, char** argv)
{
    string filename;
    int kmersize;
    int topcount;
    bool arg_ok = false;
    
    try
    {
        CmdParser cmd(argc, argv);
        
        map<string, string> params = cmd.get();
        filename = params.at(CmdParser::ARG_FILENAME);
        kmersize = stoi(params.at(CmdParser::ARG_KMERSIZE));
        topcount = stoi(params.at(CmdParser::ARG_TOPCOUNT));
        arg_ok = true;
    }
    catch(const CmdParserException& e)
    {
        cerr << e.what() << endl;
    }
    catch(const invalid_argument& e)
    {
        cerr << CmdParser::USAGE << endl;
    }
    catch(const out_of_range& e)
    {
        cerr << CmdParser::USAGE << endl;
    }

    
    if(arg_ok)
    {
        QRead qread(filename);
        qread.kmer_freq(kmersize, topcount);
    }
    
    return 0;
}