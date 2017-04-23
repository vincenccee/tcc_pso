#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>

#define TAM_POP 5
#define DIMENSION 5
#define ITERATIONS 5000
#define RUNS 5

#include "individual.hpp"
#include "population.hpp"

using namespace std;

void showPopulation();
Population *pop;

int main(int argc, char **argv) {
  pop = new Population(TAM_POP, DIMENSION, -32.0, 32.0);
  pop->initializePopulation();
  showPopulation();

  return 0;
}


void showPopulation(){
  std::vector<double> position;
  for(int i=0; i<TAM_POP; i++){
    cout << "ind " << i << endl;
    position.clear();
    position = pop->getIndividual(i)->getCurrentPosition();
    cout << "position: ";
    for(unsigned int j=0; j<position.size() ; j++){
      cout << " * " << position[j];
    }
    cout << endl << "best: ";
    position.clear();
    position = pop->getIndividual(i)->getBestPosition();
    for(unsigned int j=0; j<position.size() ; j++){
      cout << " * " << position[j];
    }
    cout << endl << "velocity: ";
    position.clear();
    position = pop->getIndividual(i)->getVelocity();
    for(unsigned int j=0; j<position.size() ; j++){
      cout << " * " << position[j];
    }
    cout << endl << " - fitness: " << pop->getIndividual(i)->getFitness() << " - Bfitness: " << pop->getIndividual(i)->getBestFitness() <<endl;
  }
}