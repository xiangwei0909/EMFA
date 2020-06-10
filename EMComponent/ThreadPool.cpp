#include "ThreadPool.h"
#include "tools.h"

#define WAIT_TIME(ms) do { std::this_thread::sleep_for(std::chrono::milliseconds(ms)); } while(0)
#define THREAD_GUARD(running, message) do { if(!running) throw std::runtime_error(message); } while(0)

using namespace component;

ThreadPool::ThreadPool() 
: running_(false), thread_num_(0)
{
}


ThreadPool::~ThreadPool()
{
}

void ThreadPool::init(size_t thread_num)
{
    auto hard_cores = getHardwareThreadNum();
    if (hard_cores != 0 && hard_cores < thread_num)
        thread_num = hard_cores;
    thread_num_ = thread_num;
    state_.resize(thread_num, State::idle);
    threads_.resize(thread_num_);
    running_ = true;
    createThread();
}

void ThreadPool::submit(TaskType task, size_t priority)
{
    assert(running_);
    THREAD_GUARD(running_, "none working thread to accept task");
    push_task(std::make_pair(priority, task));
}

void ThreadPool::run()
{
    assert(running_);
    THREAD_GUARD(running_, "none working thread to run");
    internalRun();
    while (!done()) WAIT_TIME(10);
}

void ThreadPool::interactiveRun()
{
    assert(running_);
    THREAD_GUARD(running_, "none working thread to run");
    TOOL size_t total_progress = submit_.size();
    internalRun();
    TOOL tool::BarAndPercent bar_perc;
    while (!done())
    {
        TOOL bar_perc(total_progress - remainTask(), total_progress);
        WAIT_TIME(10);
    }
}

bool ThreadPool::done() const
{
    return !hasTask() && std::all_of(state_.begin(), state_.end(), [](int val) { return val == State::idle; });
}

void ThreadPool::clear()
{
    if(running_)
        destroyThread();
}

size_t ThreadPool::getThreadNum() const
{
    return thread_num_;
}

size_t ThreadPool::getHardwareThreadNum() const
{
    return std::thread::hardware_concurrency();
}

void ThreadPool::reportInfo(Qostream & strm) const
{
    if (running_)
    {
        strm << HEADING "Thread Pool Information\n" TRAILING
            << LEVEL1 "Hardware core number:" << getHardwareThreadNum() << '\n'
            << LEVEL1 "Running thread number:" << getThreadNum() << '\n';
        strm.flush();
    }
}

bool ThreadPool::hasTask() const
{
    MutexGuard guard(mutex_);
    bool empty = tasks_.empty();
    return !empty;
}

void component::ThreadPool::push_task(PriTask task)
{
    MutexGuard guard(mutex_);
    submit_.push_back(task);
}

ThreadPool::PriTask ThreadPool::pop_task()
{
    PriTask task;
    MutexGuard guard(mutex_);
    if (!tasks_.empty())
    {
        task = tasks_.front();
        tasks_.pop_front();
    }
    return task;
}

size_t ThreadPool::remainTask() const
{
    MutexGuard guard(mutex_);
    return tasks_.size();
}

void ThreadPool::runPerThread(size_t id)
{
    while (running_)
    {
        while (hasTask())
        {
            state_[id] = State::busy;
            auto task = pop_task().second;
            if (task) task();
        }
        state_[id] = State::idle;
        std::unique_lock<std::mutex> lk(cv_mtx_);
        cv_.wait(lk);
    }
}

void ThreadPool::internalRun()
{
    auto checker = [](const PriTask& pt) { return pt.first == 0; };
    auto sorter = [](const PriTask& pt1, const PriTask& pt2) { return pt1.first > pt2.first; };
    {
        MutexGuard guard(mutex_);
        if (!std::all_of(submit_.begin(), submit_.end(), checker))
            std::sort(submit_.begin(), submit_.end(), sorter);
        tasks_.swap(submit_);
    }
    cv_.notify_all();
}

void ThreadPool::createThread()
{
    try {
        for (size_t i = 0; i < thread_num_; ++i)
            threads_[i] = std::thread(&ThreadPool::runPerThread, this, i);
    }
    catch (...)
    {
        throw std::runtime_error("error has occurred in ThreadPool::createThread()");
    }
}

void ThreadPool::destroyThread()
{
    running_ = false;
    cv_.notify_all();
    std::for_each(threads_.begin(), threads_.end(), [](std::thread& t) { t.join(); });
}
