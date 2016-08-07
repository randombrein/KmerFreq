/**
 @file
     Util.hpp
     Program that produce the most frequent DNA-k-mers of arbitrary length sorted by frequency in a given FASTQ file.
 
 @author
     Name:          Evren KANALICI
     Date:          06/08/2016
     E-Mail:        kanalici.evren@gmail.com
 */
#pragma once

#include <map>
#include <fstream>

// C-like fatal error handling macro
#define fatal_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while (0)


namespace KmerFreq {
    

class Util {

public:
    
    /**
        flips std::pair
     
        @param p    - pair to flip
        @return     - fliped std::pair
     */
    
    template<typename A, typename B>
    static std::pair<B,A> flip_pair(const std::pair<A,B> &p)
    {
        return std::pair<B,A>(p.second, p.first);
    }

    /**
        flips std::map like container to bimap (value-associated container)
     
        @param src  - std::map like container to flip
        @return     - value-associated (flipped ) std::multimap container
     */
    template<typename A, typename B, template<class,class,class...> class M, class... Args>
    static std::multimap<B,A,std::greater<B>> flip_map(const M<A,B,Args...> &src)
    {
        std::multimap<B,A,std::greater<B>> dst;
        std::transform(src.begin(), src.end(),
                       std::inserter(dst, dst.begin()),
                       flip_pair<A,B>);
        return dst;
    }

    /**
        skips lines in fstream for number of times
     
        @param fs   - file stream
        @n          - # times of skip
     */
    static void skip_line(std::ifstream& fs, int num=1)
    {
        for(int i=0; i<num; ++i)
            fs.ignore(std::numeric_limits<std::streamsize>::max(), '\n' );
    }

};
  
    
/* end namespace KmerFreq */
};
