#include "clan_particle_swarm.hpp"

using namespace std;

ClanParticleSwarm::ClanParticleSwarm(Problem *problem, int tamPopulation, int numClans){
  this->problem = problem;
  this->tamPopulation = tamPopulation;
  this->numClans = numClans;
  this->clanSize = tamPopulation/numClans;
  this->util = new Util();
}

ClanParticleSwarm::ClanParticleSwarm(){}

ClanParticleSwarm::~ClanParticleSwarm(){}

void ClanParticleSwarm::showPopulation(){
  std::vector<double> position;
  for(int i=0; i<tamPopulation; i++){
    cout << "ind " << i << "\t";
    position.clear();
    position = swarm->getIndividual(i)->getCurrentPosition();
    for(unsigned int j=0; j<position.size() ; j++){
      cout << " * " << position[j];
    }
    cout << " - fitness: " << swarm->getIndividual(i)->getFitness() << " - bestF: " << swarm->getIndividual(i)->getBestFitness() << endl;
    /* show population velocity */
    cout << "velo: " << "\t";
    for(unsigned int j=0; j<position.size() ; j++){
      cout << " * " << swarm->getIndividual(i)->getVelocity(j);
    }
    cout << endl;
  }
}

void ClanParticleSwarm::evolutionaryCicle(int iterations, int runs){
  this->iterations = iterations;
  this->runs = runs;
  this->util->openData("testdata1.txt");
  double result;

  for(int i=0; i<this->iterations; i++){
    bestPopulationFitness.push_back(0.0);
    bestIndividualFitness.push_back(0.0);
    populationDiversity.push_back(0.0);
  }
  for(int i=0; i<this->numClans; i++){
    leaders.push_back(0);
  }

  this->vMax = (problem->getUpperBound(0) - problem->getLowerBound(0)) * VMAX;

  for(int j=0; j<this->runs; j++){
    this->swarm = new Population(tamPopulation, problem->getDimension(), problem->getLowerBound(j), problem->getUpperBound(j));
    swarm->initializePopulation();
    this->m_nmdf = 0;
    evaluatePopulationFitnessFirst();
    // showPopulation();
    initializeBest();
    cout << "************** inicio ****************" << endl;
    for(int i=0; i<this->iterations; i++){
      for(int clan=0; clan<this->numClans; clan++){
        updateParticleVelocity(clan);
        updateParticlePosition(clan);
        evaluatePopulationFitness(clan);
        updateClanLeaders(clan);
      }
      leadersConference();
      // if(i%100 == 0){
      //   cout << "**** " << i << " *****" << endl;
      //   showPopulation();
      // }
      if(problem->isDynamic()){
        detectChange(i);
      }
      updateBest(i);
      updatePlot(i);
    }
    result = this->bestFitness;
    finalFitness.push_back(result);

    showPopulation();
    cout << "\nBest Finess = " << result << endl;
    problem->resetProblem();
  }
  for(int j=0; j<this->iterations; j++){
    bestPopulationFitness[j] = bestPopulationFitness[j]/this->runs;
    bestIndividualFitness[j] = bestIndividualFitness[j]/this->runs;
    populationDiversity[j] = populationDiversity[j]/this->runs;
  }
  cout << "************** Final ****************" << endl;
  cout << "Media do Fitness: " << this->util->arithmeticAverage(finalFitness) << endl;
  cout << "Desvio Padrao: " << this->util->standardDeviation(finalFitness) << endl;

  this->util->gnu_plot_convergence_best_mean(bestIndividualFitness, bestPopulationFitness, iterations, "MelhoraFitness", "melhora_fit");
  this->util->gnu_plot_convergence(populationDiversity, iterations, "pop_diversity", "fitnessdaPopulacao", "Divesidade Genotipica", 1);
  this->util->closeData();
}

void ClanParticleSwarm::updateParticleVelocity(int clan){
  std::vector<double> newVelocity;
  Individual *tmpInd;
  double cognitive, social, inertia, construtiveInertia;
  if(clan < 0){
    for(int i=0; i<numClans; i++){
      tmpInd = swarm->getIndividual(leaders[i]);
      newVelocity.clear();
      for(unsigned int j=0; j<problem->getDimension(); j++){
        if(VTYPE == 1){
          int phi = CONS1 + CONS2;
          if(phi < 4) phi = 4;
          int aux = 2 - phi - sqrt(pow(phi, 2) - 4*phi);
          if(aux < 0) aux = aux*(-1);
          construtiveInertia = 2/aux;
        } else {
          construtiveInertia = 1;
        }
        cognitive = CONS1 * this->util->fRand(0,1) * (tmpInd->getBestPosition(j) - tmpInd->getPosition(j));
        social = CONS2 * this->util->fRand(0,1) * (this->bestPosition[j] - tmpInd->getPosition(j));
        inertia = INERTIA * tmpInd->getVelocity(j);
        newVelocity.push_back(validateVelocity(construtiveInertia*(inertia + cognitive + social)));
      }
      swarm->getIndividual(leaders[i])->setVelocity(newVelocity);
    }
  } else {
    int initClan = clan*this->clanSize;
    int finalClan = (clan+1)*this->clanSize;
    for(int i=initClan; i<finalClan; i++){
      tmpInd = swarm->getIndividual(i);
      newVelocity.clear();
      for(unsigned int j=0; j<problem->getDimension(); j++){
        if(VTYPE == 1){
          int phi = CONS1 + CONS2;
          if(phi < 4) phi = 4;
          int aux = 2 - phi - sqrt(pow(phi, 2) - 4*phi);
          if(aux < 0) aux = aux*(-1);
          construtiveInertia = 2/aux;
        } else {
          construtiveInertia = 1;
        }
        cognitive = CONS1 * this->util->fRand(0,1) * (tmpInd->getBestPosition(j) - tmpInd->getPosition(j));
        social = CONS2 * this->util->fRand(0,1) * (this->bestClanPosition[clan][j] - tmpInd->getPosition(j));
        inertia = INERTIA * tmpInd->getVelocity(j);
        newVelocity.push_back(validateVelocity(construtiveInertia*(inertia + cognitive + social)));
      }
      swarm->getIndividual(i)->setVelocity(newVelocity);
    }
  }
}

double ClanParticleSwarm::validateVelocity(double velocity){
  if(velocity > vMax){
    return vMax;
  } else if(velocity < -vMax) {
    return -vMax;
  } else {
    return velocity;
  }
}

void ClanParticleSwarm::updateParticlePosition(int clan){
  std::vector<double> newPosition;
  std::vector<double> currentPosition;
  std::vector<double> velocity;
  Individual *tmpInd;

  if(clan < 0){
    for(int i=0; i<numClans; i++){
      tmpInd = swarm->getIndividual(leaders[i]);
      currentPosition = tmpInd->getCurrentPosition();
      newPosition.clear();
      for(int j=0; j<problem->getDimension(); j++){
        newPosition.push_back(currentPosition[j] + tmpInd->getVelocity(j));
      }
      swarm->getIndividual(leaders[i])->setCurrentPosition(validatePosition(newPosition));
    }
  } else {
    int initClan = clan*this->clanSize;
    int finalClan = (clan+1)*this->clanSize;
    for(int i=initClan; i<finalClan; i++){
      tmpInd = swarm->getIndividual(i);
      currentPosition = tmpInd->getCurrentPosition();
      newPosition.clear();
      for(int j=0; j<problem->getDimension(); j++){
        newPosition.push_back(currentPosition[j] + tmpInd->getVelocity(j));
      }
      swarm->getIndividual(i)->setCurrentPosition(validatePosition(newPosition));
    }
  }
}

std::vector<double> ClanParticleSwarm::validatePosition(std::vector<double> position){
  for (unsigned int i = 0; i < position.size(); i++) {
    if(position[i] > problem->getUpperBound(i)){
      position[i] = problem->getUpperBound(i);
    }else if(position[i] < problem->getLowerBound(i)){
      position[i] = problem->getLowerBound(i);
    }
  }
  return position;
}

void ClanParticleSwarm::evaluatePopulationFitnessFirst(){
  double newFitness;
  Individual *tmpInd;
  // For para ser paralizado (OPEN-MP)
  // #pragma omp parallel for
  for(int i=0; i<tamPopulation; i++){
    tmpInd = swarm->getIndividual(i);
    newFitness = problem->evaluateFitness(tmpInd->getCurrentPosition());
    tmpInd->setBestFitness(newFitness);
    tmpInd->setBestPosition(tmpInd->getCurrentPosition());
    tmpInd->setFitness(newFitness);
    swarm->updateIndividual(*tmpInd, i);
  }
}

void ClanParticleSwarm::evaluatePopulationFitness(int clan){
  double newFitness;
  Individual *tmpInd;

  if(clan < 0){
    // For para ser paralizado (OPEN-MP)
    // #pragma omp parallel for
    for(int i=0; i<numClans; i++){
      tmpInd = swarm->getIndividual(leaders[i]);
      newFitness = problem->evaluateFitness(tmpInd->getCurrentPosition());
      if(problem->fitnesIsBetter(newFitness, tmpInd->getBestFitness())) {
        tmpInd->setBestFitness(newFitness);
        tmpInd->setBestPosition(tmpInd->getCurrentPosition());
      }
      tmpInd->setFitness(newFitness);
      swarm->updateIndividual(*tmpInd, leaders[i]);
    }
  } else {
    int initClan = clan*this->clanSize;
    int finalClan = (clan+1)*this->clanSize;
    // For para ser paralizado (OPEN-MP)
    // #pragma omp parallel for
    for(int i=initClan; i<finalClan; i++){
      tmpInd = swarm->getIndividual(i);
      newFitness = problem->evaluateFitness(tmpInd->getCurrentPosition());
      if(problem->fitnesIsBetter(newFitness, tmpInd->getBestFitness())) {
        tmpInd->setBestFitness(newFitness);
        tmpInd->setBestPosition(tmpInd->getCurrentPosition());
      }
      tmpInd->setFitness(newFitness);
      swarm->updateIndividual(*tmpInd, i);
    }
  }
}

void ClanParticleSwarm::detectChange(int it){
  double newFit = problem->evaluateFitness(this->bestPosition);
  if(this->bestFitness != newFit){
    Individual *tmpInd;
    cout << "detect change!! - " << it << endl;
    showPopulation();
    this->bestFitness = newFit;
    for(int i=0; i<numClans; i++){
      tmpInd = swarm->getIndividual(leaders[i]);
      tmpInd->setBestFitness(problem->evaluateFitness(tmpInd->getFullBestPosition()));
      swarm->updateIndividual(*tmpInd, i);
    }
  }
}

void ClanParticleSwarm::initializeBest(){
  int initClan, finalClan;
  this->bestClanFitness.clear();
  this->bestClanPosition.clear();
  this->bestFitness = swarm->getIndividual(0)->getFitness();
  this->bestPosition = swarm->getIndividual(0)->getCurrentPosition();
  for(int clan=0; clan<numClans; clan++){
    initClan = clan*this->clanSize;
    finalClan = (clan+1)*this->clanSize;
    this->bestClanFitness.push_back(swarm->getIndividual(initClan)->getFitness());
    this->bestClanPosition.push_back(swarm->getIndividual(initClan)->getCurrentPosition());
    this->leaders[clan] = initClan;
    for(int i=initClan+1; i < finalClan; i++) {
      if(problem->fitnesIsBetter(swarm->getIndividual(i)->getFitness(), this->bestClanFitness[clan])){
        this->bestClanFitness[clan] = swarm->getIndividual(i)->getFitness();
        this->bestClanPosition[clan] = swarm->getIndividual(i)->getCurrentPosition();
        this->leaders[clan] = i;
        if(problem->fitnesIsBetter(this->bestClanFitness[clan], this->bestFitness)){
          this->bestFitness = this->bestClanFitness[clan];
          this->bestPosition = this->bestClanPosition[clan];
        }
      }
    }
  }
}

void ClanParticleSwarm::updateClanLeaders(int clan){
  int initClan = clan*this->clanSize;
  int finalClan = (clan+1)*this->clanSize;
  for(int i=initClan; i < finalClan; i++) {
    if(problem->fitnesIsBetter(swarm->getIndividual(i)->getFitness(), this->bestClanFitness[clan])){
      this->bestClanFitness[clan] = swarm->getIndividual(i)->getFitness();
      this->bestClanPosition[clan] = swarm->getIndividual(i)->getCurrentPosition();
      this->leaders[clan] = i;
      if(problem->fitnesIsBetter(this->bestClanFitness[clan], this->bestFitness)){
        this->bestFitness = this->bestClanFitness[clan];
        this->bestPosition = this->bestClanPosition[clan];
      }
    }
  }
}

void ClanParticleSwarm::leadersConference(){
  int clan = -1; // leader clan
  updateParticleVelocity(clan);
  updateParticlePosition(clan);
  evaluatePopulationFitness(clan);
}

void ClanParticleSwarm::updateBest(int pos){
  for(int i=0; i<numClans; i++){
    if(problem->fitnesIsBetter(swarm->getIndividual(leaders[i])->getFitness(), this->bestFitness)){
      this->bestFitness = swarm->getIndividual(leaders[i])->getFitness();
      this->bestPosition = swarm->getIndividual(leaders[i])->getCurrentPosition();
    }
  }
  // cout << "best: ";
  // for(int i=0; i < problem->getDimension(); i++) {
  //   cout << "* " << this->bestPosition[i];
  // }
  // cout << " - fitness : " << this->bestFitness << endl;
  bestIndividualFitness[pos] += this->bestFitness;
}

void ClanParticleSwarm::updatePlot(int pos){
  double totalFit = 0;
  double mediaFit = 0;

  for(int i=0; i < tamPopulation; i++) {
    totalFit += swarm->getIndividual(i)->getFitness();
  }
  mediaFit = totalFit/tamPopulation;
  bestPopulationFitness[pos] += mediaFit;
  populationDiversity[pos] += defaultGenotypicDiversityMeasure();
}

double ClanParticleSwarm::defaultGenotypicDiversityMeasure(){
  double diversity = 0;
  double aux_1 = 0;
  double aux_2 = 0;
  unsigned short int a = 0;
  unsigned short int b = 0;
  unsigned short int d = 0;
  for(a = 0; a < tamPopulation; a++)
  {
    for(b = (a+1); b < tamPopulation; b++)
    {
      aux_1 = 0;
      for(d = 0; d < problem->getDimension(); d++)
      {
        aux_1 += pow(swarm->getIndividual(a)->getCurrentPosition()[d] - swarm->getIndividual(b)->getCurrentPosition()[d], 2);
      }
      if(b == (a+1) || aux_2 > aux_1)
      {
        aux_2 = aux_1;
      }
    }
    diversity += log((double)1.0 + aux_2);
  }
  if(m_nmdf < diversity)
  {
    m_nmdf = diversity;
  }
  return diversity / m_nmdf;
}