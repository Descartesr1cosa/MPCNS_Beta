// HallConfig.h
#pragma once

#define HALL_MODE 1

#ifndef HALL_MODE
#define HALL_MODE 0 // 0: Ideal MHD, 1: Explicit Hall, 2: Implicit Hall
#endif

#if (HALL_MODE < 0) || (HALL_MODE > 2)
#error "HALL_MODE must be 0 (Ideal), 1 (Explicit Hall), or 2 (Implicit Hall)."
#endif