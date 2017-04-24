#include "rosembrock.hpp"

Rosembrock::Rosembrock(int dimension){
  this->dimension = dimension;
}

Rosembrock::Rosembrock() {}
Rosembrock::~Rosembrock() {}

double Rosembrock::getUpperBound(int pos){
  return 30.0;
}

double Rosembrock::getLowerBound(int pos){
  return -30.0;
}

double Rosembrock::evaluateFitness(std::vector<double> solution){
  double auxSum = 0;
  double aux = 0;
  for (int i = 0; i < this->dimension; i++)
  {
    aux = 100*(pow((solution[i+1] - pow(solution[1], 2)), 2));
    aux += pow(1.0 - solution[i], 2);
    auxSum += aux;
  }
  return auxSum;
}

std::string Rosembrock::getName(){
  return "Rosembrock";
}

bool Rosembrock::fitnesIsBetter(double newFit, double oldFit){
  return newFit<oldFit;
}

bool Rosembrock::isMinimization(){
  return true;
}

void Rosembrock::resetProblem() {}