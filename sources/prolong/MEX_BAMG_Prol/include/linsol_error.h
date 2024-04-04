

#pragma once

#include<string>
#include<iostream>

#define kB 1024L
#define MB 1024L*kB
#define GB 1024L*MB

#define INFO_SZ 6

#define MemUM MB

#if MemUM == kB
#define UM_label "kb"
#elif MemUM == MB
#define UM_label "Mb"
#elif MemUM == GB
#define UM_label "Gb"
#endif


class linsol_error{

   public:


   linsol_error(std::string const & __restrict__ function_name,
                std::string const & __restrict__ error_message);


   linsol_error(std::string const & __restrict__ function_name,
                std::string const & __restrict__ error_message,
                unsigned long int nbytes_requested);

};
