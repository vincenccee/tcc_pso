#include "scenario.hpp"

Scenario::Scenario(int scenario){
  switch(scenario){
    case 1:
      setScenario1();
      break;
    case 2:
      setScenario2();
      break;
    case 3:
      setScenario3();
      break;
  }
}

Scenario::Scenario(){}
Scenario::~Scenario(){}

void Scenario::setScenario1(){
  this->pfunc = 3;
  this->npeaks = 5;
  this->period = 5000;
  this->min_coord = 0.0;
  this->max_coord = 100.0;
  this->min_height = 30.0;
  this->max_height = 70.0;
  this->min_width = 0.00001;
  this->max_width = 0.2;
  this->start_height = 50.0;
  this->start_width = 0.1;
  this->move_severity = 1.0;
  this->height_severity = 7.0;
  this->width_severity = 0.01;
  this->lambda = 0.0;
}

void Scenario::setScenario2(){
  this->pfunc = 1;
  this->npeaks = 10;
  this->period = 5000;
  this->min_coord = 0.0;
  this->max_coord = 100.0;
  this->min_height = 30.0;
  this->max_height = 70.0;
  this->min_width = 1.0;
  this->max_width = 12.0;
  this->start_height = 50.0;
  this->start_width = 0.0;
  this->move_severity = 1.0;
  this->height_severity = 7.0;
  this->width_severity = 1.0;
  this->lambda = 0.5;
}

void Scenario::setScenario3(){
  this->pfunc = 1;
  this->npeaks = 5;
  this->period = 1000;
  this->min_coord = 0.0;
  this->max_coord = 100.0;
  this->min_height = 30.0;
  this->max_height = 70.0;
  this->min_width = 1.0;
  this->max_width = 12.0;
  this->start_height = 50;
  this->start_width = 0.0;
  this->move_severity = 1.0;
  this->height_severity = 1.0;
  this->width_severity = 0.5;
  this->lambda = 0.5;
}