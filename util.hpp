#ifndef _UTIL_H_
#define _UTIL_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>

#define RANDOM_DISTRIBUTION 2

using namespace std;

class Util {
  private:
    string OUTPUT_DIR;
    ofstream popdata;
    random_device rd;

  public:
    Util();
    ~Util();

    void openData(string filename);
    void closeData();
    double fRand(double fMin, double fMax);
    void gnu_plot_convergence(vector<double> mean_gen, int m_gen, string name, string title, string y_axis, double max_range);
    void gnu_plot_convergence_best_mean(vector<double> d_data1, vector<double> d_data2, int n_lines, string title, string filename);
    string space2underscore(string text);
    double standardDeviation(vector<double> data);
    double arithmeticAverage(vector<double> data);
    double euclideanDistance(vector<double> data1, vector<double> data2);
};

#endif
