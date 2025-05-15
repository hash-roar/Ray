/**
 * @file config.h 
 * @author zlf
 * @brief 
 * @version 0.1
 * @date 2025-05-15
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#pragma once

#include <thread>
namespace Ray {


    struct RayConfig{
        int seed;
        int numThreads;
        int numSamples;
        int maxDepth;
        int width;
        int height;


    };


} // namespace Ray