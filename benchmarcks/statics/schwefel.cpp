#include "schwefel.hpp"

Schwefel::Schwefel(int dimension){
  this->dimension = dimension;
}

Schwefel::Schwefel() {}
Schwefel::~Schwefel() {}

double Schwefel::getUpperBound(int pos){
  return 100.0;
}

double Schwefel::getLowerBound(int pos){
  return -100.0;
}

double Schwefel::evaluateFitness(std::vector<double> solution){
  double fitness = 0;
  double auxFit = 0;

  for (int i = 0; i < this->dimension; i++) {
    for (int j = 0; j <= i; j++) {
      auxFit += solution[j];
    }
    fitness += (auxFit * auxFit);
    auxFit = 0;
  }
  return fitness;
}

std::string Schwefel::getName(){
  return "Schwefel";
}

bool Schwefel::fitnesIsBetter(double newFit, double oldFit){
  return newFit<oldFit;
}

bool Schwefel::isMinimization(){
  return true;
}

bool Schwefel::isDynamic(){
  return false;
}

void Schwefel::resetProblem() {}