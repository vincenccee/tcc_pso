#ifndef _PARTICLE_SWARM_H_
#define _PARTICLE_SWARM_H_
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

class ParticleSwarm {
  private:
    Problem *problem;
    Population *swarm;
    std::vector<double> bestPosition;
    std::vector<double> bestPopulationFitness;
    std::vector<double> bestIndividualFitness;
    std::vector<double> populationDiversity;
    std::vector<double> finalFitness;
    std::string OUTPUT_DIR;
    std::ofstream popdata;

    double bestFitness;
    double m_nmdf;
    double velocityInertia;
    double acelerationCons1;
    double acelerationCons2;
    int runs;
    int iterations;
    int tamPopulation;
    int randonDistribution;

  public:
    ParticleSwarm(Problem *problem, int tamPopulation);
    ParticleSwarm();
    ~ParticleSwarm();

    void showPopulation();
    void evolutionaryCicle(int iterations, int runs);
    void evaluatePopulationFitness(bool first = false);
    void updateParticleVelocity();
    void updateParticlePosition();
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