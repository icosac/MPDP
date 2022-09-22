//
// Created by Paolo Bevilacqua on 05/09/17.
//

#pragma once

#include <chrono>

class TimePerf
{
public:
    TimePerf() = default;
    TimePerf(const TimePerf& tp) = default;
    virtual ~TimePerf() = default;

    void start() { startTime = std::chrono::high_resolution_clock::now(); }

    template <typename TimeGran = std::milli>
    double getTime()
    {
        auto endTime = std::chrono::high_resolution_clock::now();
        auto diff = endTime - startTime;
        return std::chrono::duration<double, TimeGran>(diff).count();
    }
    
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
};

 
