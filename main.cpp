#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

#define TAM_POP 100
#define DIMENSION 10
#define SCENARIO 2
#define ITERATIONS 2000
#define RUNS 1
#define NUM_CLANS 2

#include <string>
#include <cmath>
#include "algorithms/particle_swarm.hpp"
#include "algorithms/clan_particle_swarm.hpp"
#include "benchmarcks/problem.hpp"
#include "benchmarcks/problem_factory.hpp"

using namespace std;

string getProblem();
int getAlgorithm();

int main(int argc, char **argv) {
  string value;
  int algorithm;
  ParticleSwarm *pso;
  ClanParticleSwarm *cpso;

  value = getProblem();
  cout << "opcao escolhida foi: " << value << endl;
  ProblemFactory *factory = new ProblemFactory();
  Problem *problem = factory->get(value, DIMENSION, SCENARIO);

  algorithm = getAlgorithm();
  switch (algorithm){
    case 1:
      pso = new ParticleSwarm(problem, TAM_POP);
      pso->evolutionaryCicle(ITERATIONS, RUNS);
      break;
    case 2:
      cpso = new ClanParticleSwarm(problem, TAM_POP, NUM_CLANS);
      cpso->evolutionaryCicle(ITERATIONS, RUNS);
      break;
    default:
      cout << "Opcao Invalida" << endl;
      getAlgorithm();
  }

  return 0;
}

string getProblem(){
  int index;
  string answer;
  cout << "Select the problem:" << endl;
  cout << "1. Ackley" << endl;
  cout << "2. Rastring" << endl;
  cout << "3. Rosembrock" << endl;
  cout << "4. Griewank" << endl;
  cout << "5. Schwefel 1.2" << endl;
  cout << "6. Sphere" << endl;
  cout << "7. Moving Peaks" << endl;
  cin >> index;
  switch (index){
    case 1:
      answer = "ACKLEY";
      break;
    case 2:
      answer = "RASTRING";
      break;
    case 3:
      answer = "ROSEMBROCK";
      break;
    case 4:
      answer = "GRIEWANK";
      break;
    case 5:
      answer = "SCHWEFEL";
      break;
    case 6:
      answer = "SPHERE";
      break;
    case 7:
      answer = "MOVING_PEAKS";
      break;
    default:
      cout << "Opcao Invalida" << endl;
      getProblem();
  }
  return answer;
}

int getAlgorithm(){
  int index;
  cout << "Select the algorithm:" << endl;
  cout << "1. PSO" << endl;
  cout << "2. CPSO" << endl;
  cin >> index;
  return index;
}
