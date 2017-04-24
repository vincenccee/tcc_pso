#ifndef _ROSEMBROCK_H_
#define _ROSEMBROCK_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "../problem.hpp"

using namespace std;

class Rosembrock: public Problem {
  public:
    Rosembrock(int dimension);
    Rosembrock();
    ~Rosembrock();

    double getUpperBound(int pos);
    double getLowerBound(int pos);
    double evaluateFitness(std::vector<double> solution);
    std::string getName();
    bool fitnesIsBetter(double newFit, double oldFit);
    bool isMinimization();
    void resetProblem();
}; 

#endif