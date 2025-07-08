#ifndef _unix_cp_h
#define _unix_cp_h

#include <string>

// prototypes of functions from "unix_cp.cxx"
void copy_file(const char *file_in, const char* file_out);
void copy_dir(const char *dir_in, const char* dir_out);
void rm_whole_dir(const char *dir);
std::string copy_dir_temp(std::string dir_in, std::string basedir);
#endif
