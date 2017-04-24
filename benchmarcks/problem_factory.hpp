#ifndef PROBLEM_FACTORY_H
#define PROBLEM_FACTORY_H

#include <string>

#include "problem.hpp"
#include "statics/ackley.hpp"
#include "statics/rastring.hpp"
#include "statics/griewank.hpp"
#include "statics/rosembrock.hpp"
#include "statics/schwefel.hpp"
#include "statics/sphere.hpp"
#include "dynamics/moving_peaks.hpp"

class ProblemFactory{
public:
  Problem * get(std::string name, int dimension, int scenario = 1);
};

Problem * ProblemFactory::get(std::string name, int dimension, int scenario)
{
  if(name == "RASTRING"){
    return new Rastring(dimension);
  }
  else if(name == "ROSEMBROCK"){
    return new Rosembrock(dimension);
  }
  else if(name == "SCHWEFEL"){
    return new Schwefel(dimension);
  }
  else if(name == "GRIEWANK"){
    return new Griewank(dimension);
  }
  else if(name == "ACKLEY"){
    return new Ackley(dimension);
  }
  else if(name == "SPHERE"){
    return new Sphere(dimension);
  }
  else if(name == "MOVING_PEAKS"){
    return new MovingPeaks(dimension, scenario);
  }
  throw std::invalid_argument( "received invalid Objective Function name" );
}

#endif