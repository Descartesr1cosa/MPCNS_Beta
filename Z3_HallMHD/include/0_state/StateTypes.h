// 0_state/StateTypes.h
#pragma once
#include <cstdio>  // std::fprintf
#include <cstdlib> // std::abort

template <class FieldT>
struct Triplet
{
    FieldT *xi = nullptr;
    FieldT *eta = nullptr;
    FieldT *zeta = nullptr;

    bool all_valid() const { return xi && eta && zeta; }
    void require_all(const char *what) const
    {
        if (!all_valid())
        {
            std::fprintf(stderr, "[Triplet] %s not fully bound\n", what);
            std::abort();
        }
    }
};

// 用于缓存“字段 id 的三分量”，和 Triplet(FieldT*) 不同：这里存 int
struct IdTriplet
{
    int xi = -1;
    int eta = -1;
    int zeta = -1;

    bool all_valid() const { return (xi >= 0) && (eta >= 0) && (zeta >= 0); }

    void require_all(const char *what) const
    {
        if (!all_valid())
        {
            std::fprintf(stderr, "[IdTriplet] %s not fully bound (xi=%d eta=%d zeta=%d)\n",
                         what ? what : "(null)", xi, eta, zeta);
            std::abort();
        }
    }

    const int &at(int dir) const
    {
        switch (dir)
        {
        case 0:
            return xi;
        case 1:
            return eta;
        case 2:
            return zeta;
        default:
            std::fprintf(stderr, "[IdTriplet] invalid dir=%d (expect 0/1/2)\n", dir);
            std::abort();
        }
    }
};

// 运行时 Hall 模式（建议逐步替代 HALL_MODE 宏在主线里的散布）
enum class HallMode : int
{
    Ideal = 0,
    Explicit = 1,
    Implicit = 2
};
