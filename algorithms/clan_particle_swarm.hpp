#ifndef _CLAN_PARTICLE_SWARM_H_
#define _CLAN_PARTICLE_SWARM_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include "../population.hpp"
#include "../individual.hpp"
#include "../util.hpp"
#include "../benchmarcks/problem.hpp"

#define CONS1 2.05
#define CONS2 2.05
#define INERTIA 0.729
#define VMAX 0.1
#define VTYPE 1

#define RATE 200

class ClanParticleSwarm {
  private:
    Problem *problem;
    Population *swarm;
    Individual *testParticle;
    Util *util;
    std::vector<std::vector<double>> bestClanPosition;
    std::vector<double> bestClanFitness;
    std::vector<double> bestPosition;
    std::vector<double> bestPopulationFitness;
    std::vector<double> bestIndividualFitness;
    std::vector<double> populationDiversity;
    std::vector<double> finalFitness;
    std::vector<double> offlineError;
    std::vector<int> leaders;

    double bestFitness;
    double vMax;
    double m_nmdf;
    int numClans;
    int clanSize;
    int runs;
    int iterations;
    int tamPopulation;

  public:
    ClanParticleSwarm(Problem *problem, int tamPopulation, int numClans);
    ClanParticleSwarm();
    ~ClanParticleSwarm();

    void showPopulation();
    void contourMapData(int change);
    void evolutionaryCicle(int iterations, int runs);
    void initializeVariables();
    void initializeTestParticle();
    void evaluatePopulationFitnessFirst();
    void evaluatePopulationFitness(int clan);
    void updateParticleVelocity(int clan);
    void updateParticlePosition(int clan);
    void updateClanLeaders(int clan);
    void leadersConference();
    void detectChange(int it);
    void reevaluteBestFitness();
    double validateVelocity(double velocity);
    std::vector<double> validatePosition(std::vector<double> position);
    void initializeBest();
    void updateBest(int pos);
    void updatePlot(int pos);
    double defaultGenotypicDiversityMeasure();
};

#endif
