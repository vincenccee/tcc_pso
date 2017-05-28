#ifndef _SCHWEFEL_H_
#define _SCHWEFEL_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "../problem.hpp"

using namespace std;

class Schwefel: public Problem {
  public:
    Schwefel(int dimension);
    Schwefel();
    ~Schwefel();

    double getUpperBound(int pos);
    double getLowerBound(int pos);
    double evaluateFitness(std::vector<double> solution);
    std::string getName();
    bool fitnesIsBetter(double newFit, double oldFit);
    bool isMinimization();
    bool isDynamic();
    void resetProblem();
};

#endif