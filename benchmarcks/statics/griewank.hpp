#ifndef _GRIEWANK_H_
#define _GRIEWANK_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "../problem.hpp"

using namespace std;

class Griewank: public Problem {
  public:
    Griewank(int dimension);
    Griewank();
    ~Griewank();

    double getUpperBound(int pos);
    double getLowerBound(int pos);
    double getFitnessObjetive();
    double evaluateFitness(std::vector<double> solution);
    double evaluateFit(std::vector<double> solution);
    std::string getName();
    bool fitnesIsBetter(double newFit, double oldFit);
    bool isMinimization();
    bool isDynamic();
    void resetProblem();
};

#endif
