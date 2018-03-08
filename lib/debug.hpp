#pragma once

#include <iostream>
#include <ctime>

class timer {

public:
    timer() {
        restart();
    }
    ~timer() {
    }

    void restart() {
        start = std::clock();
    }

    double seconds_elapsed(bool restart=false) {
        double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        if (restart) this->restart();
        return duration;
    }

private:
    std::clock_t start;

};
