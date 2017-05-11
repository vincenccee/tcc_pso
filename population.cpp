#include "population.hpp"

using namespace std;

Population::Population(int tamPopulation, int dimension, double lowerBound, double upperBound){
  this->tamPopulation = tamPopulation;
  this->dimension = dimension;
  this->lowerBound = lowerBound;
  this->upperBound = upperBound;
  srand(time(NULL));
}

Population::Population(){}

Population::~Population(){}

void Population::initializePopulation(){
  this->population.clear();
  std::vector<double> position;
  std::vector<double> velocity;
  Individual *tmpInd;
  for(int i=0; i< tamPopulation; i++){
    position.clear();
    velocity.clear();
    for(int j=0; j<dimension; j++){
      position.push_back(fRand(lowerBound, upperBound));
      velocity.push_back((upperBound - lowerBound) * fRand(0,1) + lowerBound);
    }
    tmpInd = new Individual(position, velocity, this->dimension);
    this->population.push_back(*tmpInd);
  }
}

std::vector<Individual> Population::getPopulation(){
  return this->population;
}

Individual * Population::getIndividual(int pos){
  return &population[pos];
}

void Population::updateIndividual(Individual individual, int pos){
  this->population[pos] = individual;
}

int Population::getTamPopulation(){
  return this->tamPopulation;
}

double Population::fRand(double fMin, double fMax){
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}
