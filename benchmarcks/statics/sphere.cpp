#include "sphere.hpp"

Sphere::Sphere(int dimension){
  this->dimension = dimension;
}

Sphere::Sphere() {}
Sphere::~Sphere() {}

double Sphere::getUpperBound(int pos){
  return 100.0;
}

double Sphere::getLowerBound(int pos){
  return -100.0;
}

double Sphere::evaluateFitness(std::vector<double> solution){

  double top = 0;
  unsigned short int i;
  for (i = 0; i < this->dimension; i++)
  {
    top=top+solution[i]*solution[i];
  }
  return top;
}

std::string Sphere::getName(){
  return "Sphere";
}

bool Sphere::fitnesIsBetter(double newFit, double oldFit){
  return newFit<oldFit;
}

bool Sphere::isMinimization(){
  return true;
}

bool Sphere::isDynamic(){
  return false;
}

void Sphere::resetProblem() {}