#pragma once
#include <functional>
#include <limits>
class derivative{
    std::function<double(double)> _function;
    static constexpr double MIN = std::numeric_limits<double>::epsilon();
public:
    derivative(std::function<double(double)> f) : _function(f){}
    // 3점 양쪽 공식으로 구현 => 한쪽 끝으로만 미분해야 하는 경우는 곤란해지는 문제가 있긴 함.
    double operator()(double x){return (_function(x+MIN)-_function(x-MIN))/(2*MIN);}
};
