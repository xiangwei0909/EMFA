//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"

namespace component {

class ConfigInvalid : public std::invalid_argument {
public:
    explicit ConfigInvalid(const std::string& what_arg) : std::invalid_argument(what_arg) { }
    explicit ConfigInvalid(const char* what_arg) : std::invalid_argument(what_arg) { }
};

class FileError : public std::runtime_error {
public:
    explicit FileError(const std::string& what_arg) : std::runtime_error(what_arg) { }
    explicit FileError(const char* what_arg) : std::runtime_error(what_arg) { }
};

} // namespace component