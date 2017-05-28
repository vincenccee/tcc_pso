#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

class Individual {
  private:
    std::vector<double> currentPosition;
    std::vector<double> bestPosition;
    std::vector<double> velocity;
    double fitness;
    double bestFitness;

  public:
    Individual(std::vector<double> position, std::vector<double> velocity, int dimension);
    ~Individual();

    std::vector<double> getCurrentPosition();
    std::vector<double> getFullBestPosition();
    std::vector<double> getVelocityVector();
    double getBestPosition(int pos);
    double getVelocity(int pos);
    double getFitness();
    double getBestFitness();
    double getPosition(int pos);

    void setCurrentPosition(std::vector<double> currentPosition);
    void setBestPosition(std::vector<double> bestPosition);
    void setVelocity(std::vector<double> velocity);
    void setFitness(double fitness);
    void setBestFitness(double bestFitness);
};

#endif
