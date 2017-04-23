#include "individual.hpp"

using namespace std;

Individual::Individual(std::vector<double> position, int dimension){
  this->currentPosition = position;
  this->bestPosition = position;
  for(unsigned int i=0; i<dimension; i++){
    this->velocity.push_back(0.0);
  }
}

Individual::~Individual(){}

std::vector<double> Individual::getCurrentPosition(){
  return this->currentPosition;
}

std::vector<double> Individual::getBestPosition(){
  return this->bestPosition;
}

std::vector<double> Individual::getVelocity(){
  return this->velocity;
}

double Individual::getFitness(){
  return this->fitness;
}

double Individual::getBestFitness(){
  return this->bestFitness;
}

double Individual::getPosition(int pos){
  return this->currentPosition[pos];
}

void Individual::setCurrentPosition(std::vector<double> currentPosition){
  this->currentPosition = currentPosition;
}

void Individual::setBestPosition(std::vector<double> bestPosition){
  this->bestPosition = bestPosition;
}

void Individual::setVelocity(std::vector<double> velocity){
  this->velocity = velocity;
}

void Individual::setFitness(double fitness){
  this->fitness = fitness;
}
