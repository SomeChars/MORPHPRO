#ifndef COMMON_H_
#define COMMON_H_

#define NDEBUG

#include<assert.h>
#include<stdio.h>
#include<stdlib.h>

#define  INF  1e+10
#define  EPS  1e-8
#define  M_PI 3.14159265358979

#define sqr(X) ((X) * (X))

extern void init_log();
extern void write_to_log(char *message);

#endif // COMMON_H_