#pragma once

#include <map>
#include <fstream>

#define handle_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while (0)


namespace KmerFreq {
    

class Util {

public:

    template<typename A, typename B>
    static std::pair<B,A> flip_pair(const std::pair<A,B> &p)
    {
        return std::pair<B,A>(p.second, p.first);
    }

    template<typename A, typename B, template<class,class,class...> class M, class... Args>
    static std::multimap<B,A,std::greater<B>> flip_map(const M<A,B,Args...> &src)
    {
        std::multimap<B,A,std::greater<B>> dst;
        std::transform(src.begin(), src.end(),
                       std::inserter(dst, dst.begin()),
                       flip_pair<A,B>);
        return dst;
    }

    static void skip_line(std::ifstream& fs, int num=1)
    {
        for(int i=0; i<num; ++i)
            fs.ignore(std::numeric_limits<std::streamsize>::max(), '\n' );
    }

};
  
    
/* end namespace KmerFreq */
};
