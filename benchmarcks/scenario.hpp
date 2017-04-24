#ifndef _SCENARIO_H_
#define _SCENARIO_H_
#include <string>

class Scenario {
  public:
    Scenario(int scenario);
    Scenario();
    ~Scenario();
    void setScenario1();
    void setScenario2();
    void setScenario3();

    int pfunc;
    int npeaks;
    int period;
    double min_coord;
    double max_coord;
    double min_height;
    double max_height;
    double min_width;
    double max_width;
    double start_height;
    double start_width;
    double move_severity;
    double height_severity;
    double width_severity;
    double lambda;
};
#endif
