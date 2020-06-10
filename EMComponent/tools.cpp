#include "tools.h"
#include <filesystem>

namespace fs = std::experimental::filesystem;
namespace ct = std::chrono;

static  std::chrono::time_point<std::chrono::system_clock>  start;
static  std::chrono::time_point<std::chrono::system_clock>  prev;

//  format time
//  @return (xxh:xxm:xxs)
Qstring tool::formatTime(long long ms)
{
    auto hours = ms / 3600000;
    ms -= hours * 3600000;
    auto minutes = ms / 60000;
    ms -= minutes * 60000;
    auto seconds = ms / 1000;
    return '(' + std::to_string(hours) + "h:" + std::to_string(minutes) + "m:" + std::to_string(seconds) + "s)";
}

Qstring tool::formatMemory(size_t mb)
{
    auto tb = mb / 1048576;
    mb -= tb * 1048576;
    auto gb = mb / 1024;
    mb -= gb * 1024;
    if (gb == 0 && tb == 0)
        return "";
    return '(' + (tb != 0 ? (std::to_string(tb) + "TB:") : "") +
            (gb != 0 ? (std::to_string(gb) + "GB:") : "") + std::to_string(mb) + "MB)";
}

void tool::initTime()
{
    start = std::chrono::system_clock::now();
    prev = start;
}

//  time cost from previous process end to current process end
//  description : cost time(ms)
//  return time cost
long long tool::elapsedTime(Qostream& strm,const Qstring& description)
{
    auto old_state = strm.flags();
    strm << std::left;
    auto current = ct::system_clock::now();
    auto period = ct::duration_cast<ct::milliseconds>(current - prev).count();
    strm << std::setw(30) << description+':' << std::right 
        << std::setw(12) << period << " ms" << std::setw(15) << formatTime(period) << '\n';
    strm.flags(old_state);
    prev = current;
    return period;
}

//  total runtime cost
long long tool::totalTime(Qostream & strm)
{
    auto old_state = strm.flags();
    strm << std::left;
    auto end = ct::system_clock::now();
    auto period = ct::duration_cast<ct::milliseconds>(end - start).count();
    strm <<std::setw(30) << "Total time:" << std::right
        << std::setw(12) << period << " ms" << std::setw(15) << formatTime(period) << '\n';
    strm.flags(old_state);
    return period;
}

//  creat folder(absolute path)
//  @return dir/foldername
Qstring tool::creatFolder(const Qstring & dir, const Qstring & foldername)
{
    auto directory = dir + '/' + foldername;
    fs::create_directory(fs::v1::path(directory));
    return directory;
}

tool::Percent::Percent() 
: old_(0)
{
    Qcout << " 0%" << std::flush;
}

tool::Percent::~Percent()
{
}

void tool::Percent::operator()(size_t cur, size_t total)
{
    size_t num = 100 * cur / total;
    if (num != old_)
    {
        old_ = num;
        Qcout << "\b\b\b" << (num < 10 ? " " : "") << num << '%' << std::flush;
    }
}

tool::ProgressBar::ProgressBar() 
: old_(0), str_("[....................]")
{
    Qcout << str_ << std::flush;
}

tool::ProgressBar::~ProgressBar()
{
}

void tool::ProgressBar::operator()(size_t cur, size_t total)
{
    size_t num = 20 * cur / total;
    if (num != old_)
    {
        if (old_ + 1 != num)
        {
            for (size_t i = 1; i < num; ++i)
                str_[i] = '#';
        }
        old_ = num;
        str_[num] = '#';
        Qcout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << str_ << std::flush;
    }
}

tool::BarAndPercent::BarAndPercent() 
: old_per_(0), old_bar_(0), str_("[....................]")
{
    Qcout << str_ << "  0%" << std::flush;
}

tool::BarAndPercent::~BarAndPercent()
{
}

void tool::BarAndPercent::operator()(size_t cur, size_t total)
{
    size_t num_per = 100 * cur / total;
    size_t num_bar = 20 * cur / total;
    if (num_per != old_per_ || num_bar != old_bar_)
    {
        old_per_ = num_per;
        Qcout << "\b\b\b\b";
        if (num_bar != old_bar_)
        {
            if (old_bar_ + 1 != num_bar)
                for (size_t i = 1; i < num_bar; ++i)
                    str_[i] = '#';
            old_bar_ = num_bar;
            str_[num_bar] = '#';
            Qcout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << str_;
        }
        Qcout << (num_per < 10 ? "  " : " ") << num_per << '%' << std::flush;
    }
}
