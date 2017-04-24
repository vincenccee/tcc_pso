#ifndef _PROBLEM_H_
#define _PROBLEM_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

class Problem {
  public:
    virtual double getUpperBound(int pos) = 0;
    virtual double getLowerBound(int pos) = 0;
    virtual double evaluateFitness(std::vector<double> solution) = 0;
    virtual std::string getName() = 0;
    virtual bool fitnesIsBetter(double newFit, double oldFit) = 0;
    virtual bool isMinimization() = 0;
    virtual void resetProblem() = 0;
    
    std::vector<double> validatePosition(std::vector<double> position);
    double getDimension();

  protected:
  	int dimension;
};

#endif
