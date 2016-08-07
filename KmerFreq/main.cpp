/**
 @file
    main.cpp
    Program that produce the most frequent DNA-k-mers of arbitrary length sorted by frequency in a given FASTQ file.

 @author
    Name:          Evren KANALICI
    Date:          06/08/2016
    E-Mail:        kanalici.evren@gmail.com
*/
#include <iostream>
#include <string>
#include <exception>

#include "CmdParser.hpp"
#include "QRead.hpp"

using namespace std;
using namespace KmerFreq;

int main(int argc, char** argv)
{
    // command arguments
    string filename;
    int kmersize;
    int topcount;
    
    // args ok to run program
    bool arg_ok = false;
    
    try
    {
        // get command line arguments
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

    // run alrorith with passed arguments
    if(arg_ok)
    {
        QRead qread(filename);
        qread.kmer_freq(kmersize, topcount);
    }
    
    return 0;
}