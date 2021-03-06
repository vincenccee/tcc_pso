#include "individual.hpp"

using namespace std;

Individual::Individual(std::vector<double> position, std::vector<double> velocity, int dimension){
  this->currentPosition = position;
  this->bestPosition = position;
  this->velocity = velocity;
  this->fitness = 0.0;
}

Individual::~Individual(){}

std::vector<double> Individual::getCurrentPosition(){
  return this->currentPosition;
}

std::vector<double> Individual::getVelocityVector(){
  return this->velocity;
}

std::vector<double> Individual::getFullBestPosition(){
  return this->bestPosition;
}

double Individual::getBestPosition(int pos){
  return this->bestPosition[pos];
}

double Individual::getVelocity(int pos){
  return this->velocity[pos];
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

void Individual::setBestFitness(double bestFitness){
  this->bestFitness = bestFitness;
}
