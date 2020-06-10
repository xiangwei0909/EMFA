#include "Test.h"

#define TEST_NAME(name) do { Qcout << HEADING "[Test Case]: " << name << "\n" TRAILING << std::flush; } while(0)

void test::Test()
{
    test::qmat::OperatorAdd();
    test::threadpool::TestThreadPool();
}

void test::qmat::OperatorAdd()
{
    TEST_NAME("test::qmat::OperatorAdd()");
    wt::QMat<size_t> m1, m2;
    m1.set_size(4, 4);
    m2.set_size(4, 4);
    for (size_t i = 0; i < m1.size(); ++i)
    {
        m1(i) = i;
        m2(i) = 10 * i;
    }
    m1 += m2;
    Qcout << "size: " << m1.size() << std::endl;
    for (size_t i = 0; i < m1.size(); ++i)
        Qcout << i << ": " << m1(i) << std::endl;
}

void test::threadpool::TestThreadPool()
{
    TEST_NAME("test::threadpool::TestThreadPool()");
    std::vector<int> nums(2, 0);
    std::function<void()> task1 =
        []() { std::this_thread::sleep_for(std::chrono::seconds(1)); Qcout << "task1: " << std::this_thread::get_id() << std::endl; };
    std::function<void()> task2 =
        []() { std::this_thread::sleep_for(std::chrono::seconds(2)); Qcout << "task2: " << std::this_thread::get_id() << std::endl; };
    std::function<void()> task3 =
        []() { std::this_thread::sleep_for(std::chrono::seconds(3)); Qcout << "task3: " << std::this_thread::get_id() << std::endl; };
    std::function<void()> task4 =
        []() { std::this_thread::sleep_for(std::chrono::seconds(4)); Qcout << "task4: " << std::this_thread::get_id() << std::endl; };
    std::function<void()> task5 =
        [&nums]() { ++nums[0]; Qcout << "task5: " << std::this_thread::get_id() << " number: " << nums[0] << std::endl; };
    std::function<void()> task6 =
        [nums]() mutable { ++nums[1]; Qcout << "task5: " << std::this_thread::get_id() << " number: " << nums[1] << std::endl; };

    component::ThreadPool pool;
    pool.init(4);
    Qcout << "thread number: " << pool.getThreadNum() << std::endl;
    Qcout << "hardware cores: " << pool.getHardwareThreadNum() << std::endl;
    pool.submit(task1);
    pool.submit(task2);
    pool.submit(task3);
    pool.submit(task4);
    pool.submit(task5);
    pool.submit(task6);
    pool.run();
    Qcout << "number1: " << nums[0] << ", number 2: " << nums[1] << std::endl;
    pool.clear();
}
