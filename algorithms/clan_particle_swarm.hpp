#ifndef _CLAN_PARTICLE_SWARM_H_
#define _CLAN_PARTICLE_SWARM_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "../population.hpp"
#include "../individual.hpp"
#include "../benchmarcks/problem.hpp"

#define CONS1 2.05
#define CONS2 2.05
#define INERTIA 0.729
#define VMAX 0.1
#define VTYPE 1

class ClanParticleSwarm {
  private:
    Problem *problem;
    Population *swarm;
    std::vector<std::vector<double>> bestClanPosition;
    std::vector<double> bestClanFitness;
    std::vector<double> bestPosition;
    std::vector<double> bestPopulationFitness;
    std::vector<double> bestIndividualFitness;
    std::vector<double> populationDiversity;
    std::vector<double> finalFitness;
    std::vector<int> leaders;
    std::string OUTPUT_DIR;
    std::ofstream popdata;

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
    void evolutionaryCicle(int iterations, int runs);
    void evaluatePopulationFitnessFirst();
    void evaluatePopulationFitness(int clan);
    void updateParticleVelocity(int clan);
    void updateParticlePosition(int clan);
    void updateClanLeaders(int clan);
    void leadersConference();
    double validateVelocity(double velocity);
    std::vector<double> validatePosition(std::vector<double> position);

    void initializeBest();
    void updateBest(int pos);
    double fRand(double fMin, double fMax);
    void updatePlot(int pos);
    void gnu_plot_convergence(std::vector<double> mean_gen, int m_gen, std::string name, std::string title, std::string y_axis, double max_range);
    void gnu_plot_convergence_best_mean(std::vector<double> d_data1, std::vector<double> d_data2, int n_lines, std::string title, std::string filename);
    std::string space2underscore(std::string text);
    double defaultGenotypicDiversityMeasure();
    double standardDeviation(std::vector<double> data);
    double arithmeticAverage(std::vector<double> data);
};

#endif