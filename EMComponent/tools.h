//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#ifndef _TOOLS_H_
#define _TOOLS_H_
#include "Custom.h"

#define TOOL
#define FORMAT_MEMORY(mem) (" " + tool::formatMemory(mem))
namespace tool {

static  Qstring formatTime(long long ms);
Qstring     formatMemory(size_t mb);

void        initTime();
long long   elapsedTime(Qostream& strm, const Qstring& description = "previous time cost");
long long   totalTime(Qostream& strm);

Qstring     creatFolder(const Qstring& dir,const Qstring& foldername);


class Percent {
public:
    Percent();
    ~Percent();
    void operator()(size_t cur, size_t total);
private:
    size_t  old_;
};

class ProgressBar {
public:
    ProgressBar();
    ~ProgressBar();
    void operator()(size_t cur, size_t total);
private:
    size_t  old_;
    Qstring str_;
};

class BarAndPercent {
public:
    BarAndPercent();
    ~BarAndPercent();
    void operator()(size_t cur, size_t total);
private:
    size_t  old_per_;
    size_t  old_bar_;
    Qstring str_;
};


} // namespace tool

#endif