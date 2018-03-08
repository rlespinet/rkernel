#pragma once

static inline int ipow(int base, int exp) {
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}

static inline int ilog(int n, int base) {

    int result = 0;
    while (n /= base) {
        result++;
    }

    return result;
}
