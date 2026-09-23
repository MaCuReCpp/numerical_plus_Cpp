#ifndef GUARD_ROOT_H
#define GUARD_ROOT_H

#include <functional>

double Bisection(std::function<double(double)>, double, double);
double Secant(std::function<double(double)>, double, double);
double Falseposition(std::function<double(double)>, double, double);

#endif
