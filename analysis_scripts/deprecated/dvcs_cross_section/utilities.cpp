// utilities.cpp
#include "utilities.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>

void ensure_directory_exists(const std::string &path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        if (mkdir(path.c_str(), 0777) != 0) {
            std::cerr << "Error: Could not create directory " << path << std::endl;
        }
    }
}