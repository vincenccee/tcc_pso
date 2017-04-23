#ifndef _POPULATION_H_
#define _POPULATION_H_
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "individual.hpp"

class Population {
  private:
    std::vector<Individual> population;
    std::vector<double> globalBest;
    double upperBound;
    double lowerBound;
    int dimension;
    int tamPopulation;

  public:
    Population(int tamPopulation, int dimension, double lowerBound, double upperBound);
    Population();
    ~Population();

    void initializePopulation();
    std::vector<double> randonPosition();
    std::vector<Individual> getPopulation();
    Individual * getIndividual(int pos);
    void updateIndividual(Individual individual, int pos);

    int getTamPopulation();
    std::vector<double> getGlobalBest();
    void setGlobalBest(std::vector<double> globalBest);

    double fRand(double fMin, double fMax);
};

#endif