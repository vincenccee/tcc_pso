#ifndef _SPHERE_H_
#define _SPHERE_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "../problem.hpp"

using namespace std;

class Sphere: public Problem {
  public:
    Sphere(int dimension);
    Sphere();
    ~Sphere();

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
