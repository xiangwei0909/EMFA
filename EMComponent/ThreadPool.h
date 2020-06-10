//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

// Thread Pool Class	
// Author:  ZhangQiang
// Date:    2017-05-12 15:09:37

#pragma once
#include "Custom.h"

namespace component {

class ThreadPool {
    using TaskType = std::function<void()>;
    using PriTask = std::pair<size_t, TaskType>;
    using MutexGuard = std::lock_guard<std::mutex>;
    enum State {idle, busy};
public:
    ThreadPool();
    ~ThreadPool();
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    void init(size_t thread_num);
    void submit(TaskType task, size_t priority = 0);    // the large value represents the high priority
    void run();
    void interactiveRun();
    void clear();
    size_t getThreadNum() const;
    size_t getHardwareThreadNum() const;
    void reportInfo(Qostream& strm) const;
private:
    bool done() const;
    bool hasTask() const;
    void push_task(PriTask task);
    PriTask pop_task();
    size_t remainTask() const;
    void runPerThread(size_t id);
    void internalRun();
    void createThread();
    void destroyThread();
private:
    bool running_;
    size_t thread_num_;
    mutable std::mutex mutex_, cv_mtx_;
    mutable std::condition_variable cv_;
    std::deque<PriTask> tasks_, submit_;
    std::vector<State> state_;
    std::vector<std::thread> threads_;
};

} // namespace component
