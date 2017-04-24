#ifndef _MOVING_PEAKS_H_
#define _MOVING_PEAKS_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

// Source: https://github.com/DEAP/deap/blob/master/deap/benchmarks/movingpeaks.py
#include "../problem.hpp"
#include "../scenario.hpp"

using namespace std;

class MovingPeaks: public Problem {
  public:
    MovingPeaks(int dimension, int scenario = 1);
    MovingPeaks();
    ~MovingPeaks();

    double getUpperBound(int pos);
    double getLowerBound(int pos);
    double evaluateFitness(std::vector<double> solution);
    std::string getName();
    bool fitnesIsBetter(double newFit, double oldFit);
    bool isMinimization();
    void resetProblem();

  private:
    void showPeakes();
    void initializePeaks();
    double cone(std::vector<double> individual, std::vector<double> position, double height, double width);
    double sphere(std::vector<double> individual, std::vector<double> position, double height, double width);
    double function1(std::vector<double> individual, std::vector<double> position, double height, double width);
    double callFunction(std::vector<double> individual, int pos);
    double maxValue(std::vector<double> values);
    void changePeaks();

    double fRand(double fMin, double fMax);

    int count;
    Scenario * scen;
    std::vector<double> peaks_height;
    std::vector<double> peaks_width;
    std::vector<std::vector<double>> last_change_vector;
    std::vector<std::vector<double>> peaks_position;
};

#endif

// The Moving Peaks Benchmark is a fitness function changing over time. It
// consists of a number of peaks, changing in height, width and location. The
// peaks function is given by *pfunc*, wich is either a function object or a
// list of function objects (the default is :func:`function1`). The number of
// peaks is determined by *npeaks* (which defaults to 5). This parameter can
// be either a integer or a sequence. If it is set to an integer the number
// of peaks won't change, while if set to a sequence of 3 elements, the
// number of peaks will fluctuate between the first and third element of that
// sequence, the second element is the inital number of peaks. When
// fluctuating the number of peaks, the parameter *number_severity* must be
// included, it represents the number of peak fraction that is allowed to
// change. The dimensionality of the search domain is *dim*. A basis function
// *bfunc* can also be given to act as static landscape (the default is no
// basis function). The argument *random* serves to grant an independent
// random number generator to the moving peaks so that the evolution is not
// influenced by number drawn by this object (the default uses random
// functions from the Python module :mod:`random`). Various other keyword
// parameters listed in the table below are required to setup the benchmark,
// default parameters are based on scenario 1 of this benchmark.

// =================== ============================= =================== =================== ======================================================================================================================
// Parameter           :data:`SCENARIO_1` (Default)  :data:`SCENARIO_2`  :data:`SCENARIO_3`    Details
// =================== ============================= =================== =================== ======================================================================================================================
// ``pfunc``           :func:`function1`             :func:`cone`        :func:`cone`        The peak function or a list of peak function.
// ``npeaks``          5                             10                  50                  Number of peaks. If an integer, the number of peaks won't change, if a sequence it will fluctuate [min, current, max].
// ``bfunc``           :obj:`None`                   :obj:`None`         ``lambda x: 10``    Basis static function.
// ``min_coord``       0.0                           0.0                 0.0                 Minimum coordinate for the centre of the peaks.
// ``max_coord``       100.0                         100.0               100.0               Maximum coordinate for the centre of the peaks.
// ``min_height``      30.0                          30.0                30.0                Minimum height of the peaks.
// ``max_height``      70.0                          70.0                70.0                Maximum height of the peaks.
// ``uniform_height``  50.0                          50.0                0                   Starting height for all peaks, if ``uniform_height <= 0`` the initial height is set randomly for each peak.
// ``min_width``       0.0001                        1.0                 1.0                 Minimum width of the peaks.
// ``max_width``       0.2                           12.0                12.0                Maximum width of the peaks
// ``uniform_width``   0.1                           0                   0                   Starting width for all peaks, if ``uniform_width <= 0`` the initial width is set randomly for each peak.
// ``lambda_``         0.0                           0.5                 0.5                 Correlation between changes.
// ``move_severity``   1.0                           1.5                 1.0                 The distance a single peak moves when peaks change.
// ``height_severity`` 7.0                           7.0                 1.0                 The standard deviation of the change made to the height of a peak when peaks change.
// ``width_severity``  0.01                          1.0                 0.5                 The standard deviation of the change made to the width of a peak when peaks change.
// ``period``          5000                          5000                1000                Period between two changes.
// =================== ============================= =================== =================== ======================================================================================================================
