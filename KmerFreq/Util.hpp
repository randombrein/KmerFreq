#pragma once

#include <map>
#include <fstream>

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

    static size_t count_line(std::istream &is)
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

    static  void skip_line(std::ifstream& fs, int num=1)
    {
        for(int i=0; i<num; ++i)
            fs.ignore(std::numeric_limits<std::streamsize>::max(), '\n' );
    }

};
