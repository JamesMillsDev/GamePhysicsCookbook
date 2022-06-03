#pragma once

#define ABSOLUTE(x, y) (fabsf(x-y) <= FLT_EPSILON)
#define RELATIVE(x, y) (fabsf(x-y) <= FLT_EPSILON * fmaxf(fabsf(x), fabsf(y)))

#define RAD2DEG(x) ((x) * 57.295754f)
#define DEG2RAD(x) ((x) * 0.0174533f)

#define CMP(x, y) (fabsf(x-y) <= FLT_EPSILON * fmaxf(1.0f, fmaxf(fabsf(x), fabsf(y))))