#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iostream>
#include <thread>

struct Timer
{
    std::chrono::_V2::system_clock::time_point start, end;
    std::chrono::duration<float> duration;

    Timer() { start = std::chrono::high_resolution_clock::now(); }

    ~Timer()
    {
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;

        float ms = duration.count() * 1000.0f;

        std::clog << "Duration: " << ms << "ms.\n";
    }
};

#endif