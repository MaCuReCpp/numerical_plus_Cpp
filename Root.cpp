#include "Root.hpp"
#include <exception>
#include <cmath>

#define MAX 10000
#define MIN 1.0e-10

double Bisection(std::function<double(double)> f, double a, double b){
    if(f(a)*f(b) > 0)
        throw std::domain_error("f(a) * f(b) is plus\n");
    double p = (a+b)/2;
    for(int i = 0; i < MAX && std::abs(b - a) > MIN; ++i){
        if(f(a) * f(p) < 0){
            b = p;
            p = (a + b) / 2;
        }
        else{
            a = p;
            p = (a + b) / 2;
        }
    }
    return p;
}

double Secant(std::function<double(double)> f, double a, double b){
    if(f(a) * f(b) > 0)
        throw std::domain_error("f(a) * f(b) > 0\n");
    if(std::abs(f(a)) >= std::abs(f(b))){
        double p1 = a;
        double p2 = b;
        for(int i = 0; i < MAX && std::abs(f(p2) - f(p1)) > MIN; ++i){
            double temp = p2;
            p2 = p2 - f(p2) * (p2 - p1) / (f(p2) - f(p1));
            p1 = temp;
        }
        return p2;
    }
    else{
        double p1 = b;
        double p2 = a;
        for(int i = 0; i < MAX && std::abs(f(p1) - f(p2)) > MIN; ++i){
            double temp = p2;
            p2 = p2 - f(p2) * (p2 - p1) / (f(p2) - f(p1));
            p1 = temp;
        }
        return p2;
    }
}

double FalsePosition(std::function<double(double)> f, double a, double b){
    if(f(a) * f(b) > 0)
        throw std::domain_error("f(a) * f(b) > 0\n");
    double p = a - (f(a) * (b - a)) / (f(b) - f(a));
    for(int i = 0; i < MAX && std::abs(f(a) - f(b)) > MIN; ++i){
        if(f(a) * f(p) < 0)
            b = p;
        else
            a = p;
        p = a - (f(a) * (b - a)) / (f(b) - f(a));
    }   
    return p;
}
