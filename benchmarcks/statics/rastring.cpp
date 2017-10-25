#include "rastring.hpp"

Rastring::Rastring(int dimension){
  this->dimension = dimension;
}

Rastring::Rastring() {}
Rastring::~Rastring() {}

double Rastring::getUpperBound(int pos){
  return 5.12;
}

double Rastring::getLowerBound(int pos){
  return -5.12;
}

double Rastring::getFitnessObjetive(){
  return 0.0;
}

double Rastring::evaluateFitness(std::vector<double> solution){
  double obj = 0;
  for(int j = 0 ; j < this->dimension; j++) {
    obj += (pow(solution[j], 2) - 10 * cos(2 * M_PI*solution[j]) + 10);
  }
  return obj;
}

double Rastring::evaluateFit(std::vector<double> solution){
  double obj = 0;
  for(int j = 0 ; j < this->dimension; j++) {
    obj += (pow(solution[j], 2) - 10 * cos(2 * M_PI*solution[j]) + 10);
  }
  return obj;
}

std::string Rastring::getName(){
  return "Rastring";
}

bool Rastring::fitnesIsBetter(double newFit, double oldFit){
  return newFit<oldFit;
}

bool Rastring::isMinimization(){
  return true;
}

bool Rastring::isDynamic(){
  return false;
}

void Rastring::resetProblem() {}
