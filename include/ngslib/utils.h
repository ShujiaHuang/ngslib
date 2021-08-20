// Author: Shujia Huang
// Date: 2021-08-18
#ifndef __INCLUDE_NGSLIB_UTILS_H__
#define __INCLUDE_NGSLIB_UTILS_H__

#include <string>
#include <sstream>
#include <unistd.h>


namespace ngslib {

    // Template function can only be defined in C++ header file
    template<typename T>
    std::string tostring(T d) {
        std::stringstream ss;
        ss << d;
        return ss.str();
    }

    bool is_readable(const char *name);
}  // namespace ngslib

#endif  // #ifndef __INCLUDE_NGSLIB_UTILS_H__
