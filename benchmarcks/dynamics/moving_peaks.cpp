#include "moving_peaks.hpp"

MovingPeaks::MovingPeaks(int dimension, int scenario){
  this->dimension = dimension;
  this->scen = new Scenario(scenario);
  initializePeaks();
}

void MovingPeaks::initializePeaks(){
  this->count = 0;
  
  // set initial height for all peaks
  if(scen->start_height != 0){
    for(int i=0; i<scen->npeaks; i++)
      peaks_height.push_back(scen->start_height);
  }else{
    for(int i=0; i<scen->npeaks; i++)
      peaks_height.push_back(fRand(scen->min_height, scen->max_height));
  }

  // set initial width for all peaks 
  if(scen->start_width != 0){
    for(int i=0; i<scen->npeaks; i++)
      peaks_width.push_back(scen->start_width);
  }else{
    for(int i=0; i<scen->npeaks; i++)
      peaks_width.push_back(fRand(scen->min_width, scen->max_width));
  }

  // set initial position and last_change_vector for all peaks
  std::vector<double> tmp_pos, last_pos;
  for(int i=0; i<scen->npeaks; i++){
    tmp_pos.clear();
    last_pos.clear();
    for(int j=0; j<this->dimension; j++){
      tmp_pos.push_back(fRand(scen->min_coord, scen->max_coord));
      last_pos.push_back(fRand(0,1) - 0.5);
    }
    peaks_position.push_back(tmp_pos);
    last_change_vector.push_back(last_pos);
  }
  showPeakes();
}

MovingPeaks::MovingPeaks(){}
MovingPeaks::~MovingPeaks(){}

void MovingPeaks::showPeakes() {
  cout << "************" << endl;
  for (int i = 0; i < scen->npeaks; ++i){
    cout << "Peak " << i+1 << ": " << endl;
    cout << "Height: " << peaks_height[i] << endl;
    cout << "Width: " << peaks_width[i] << endl;
    cout << "Position: ";
    for (int j = 0; j < this->dimension; ++j){
      cout << " - " << peaks_position[i][j];
    }
    cout << endl;
  }
  cout << "************" << endl;
}

double MovingPeaks::getUpperBound(int pos){
	return scen->max_coord;
}

double MovingPeaks::getLowerBound(int pos){
	return scen->min_coord;
}

double MovingPeaks::evaluateFitness(std::vector<double> solution){
  std::vector<double> possible_values;
  this->count++;
  for(int i=0; i<scen->npeaks; i++)
    possible_values.push_back(callFunction(solution, i));

  if(this->count % scen->period == 0){
    cout << "evaluations: " << this->count << endl;
    changePeaks();
  }

  return maxValue(possible_values);
}

std::string MovingPeaks::getName(){
	return "MOVING_PEAKS";
}

bool MovingPeaks::fitnesIsBetter(double newFit, double oldFit){
	return newFit>oldFit;
}

bool MovingPeaks::isMinimization(){
	return false;
}

double MovingPeaks::cone(std::vector<double> individual, std::vector<double> position, double height, double width){
  double value = 0.0;
  for(int i=0; i<this->dimension; i++)
    value += pow(individual[i] - position[i], 2);
  return height - width * sqrt(value);
}

double MovingPeaks::sphere(std::vector<double> individual, std::vector<double> position, double height, double width){
  double value = 0.0;
  for(int i=0; i<this->dimension; i++)
    value += pow(individual[i] - position[i], 2);
  return height  * value;
}

double MovingPeaks::function1(std::vector<double> individual, std::vector<double> position, double height, double width){
  double value = 0.0;
  for(int i=0; i<this->dimension; i++)
    value += pow(individual[i] - position[i], 2);
  return height / (1 + width * value);
}

double MovingPeaks::callFunction(std::vector<double> individual, int pos){
  double result;
  switch(scen->pfunc){
    case 1:
      result = cone(individual, peaks_position[pos], peaks_height[pos], peaks_width[pos]);
      break;
    case 2:
      result = sphere(individual, peaks_position[pos], peaks_height[pos], peaks_width[pos]);
      break;
    case 3:
      result = function1(individual, peaks_position[pos], peaks_height[pos], peaks_width[pos]);
      break;
  }
  return result;
}

void MovingPeaks::changePeaks(){
  cout << "Mudei!" << endl;
  double shift_length, new_coord, change, new_value;
  std::vector<double> shift;
  std::vector<double> new_position;
  std::vector<double> final_shift;

  for(int j=0; j<scen->npeaks; j++){
    // changing peaks position
    shift_length = new_coord = 0.0;
    shift.clear();
    new_position.clear();
    final_shift.clear();

    for(int i=0; i<this->dimension; i++){
      shift.push_back(fRand(0,1) - 0.5);
      shift_length += pow(shift[i],2);
    }
    if(shift_length > 0)
      shift_length = scen->move_severity / sqrt(shift_length);
    else
      shift_length = 0.0;
    for(int i=0; i<this->dimension; i++)
      shift[i] = shift_length * (1.0 - scen->lambda) * shift[i] + scen->lambda * last_change_vector[j][i];
    shift_length = 0.0;
    for(int i=0; i<this->dimension; i++)
      shift_length += pow(shift[i],2);
    if(shift_length > 0)
      shift_length = scen->move_severity / sqrt(shift_length);
    else
      shift_length = 0.0;
    for(int i=0; i<this->dimension; i++)
      shift[i] = shift[i] * shift_length;

    for(int i=0; i<this->dimension; i++){
      new_coord = peaks_position[j][i] + shift[i];
      if(new_coord < scen->min_coord){
        new_position.push_back(2.0 * scen->min_coord - peaks_position[j][i] - shift[i]);
        final_shift.push_back(-1.0 * shift[i]);
      }else if(new_coord > scen->max_coord){
        new_position.push_back(2.0 * scen->max_coord - peaks_position[j][i] - shift[i]);
        final_shift.push_back(-1.0 * shift[i]);
      }else{
        new_position.push_back(new_coord);
        final_shift.push_back(shift[i]);
      }
    }

    peaks_position[j] = new_position;
    last_change_vector[j] = final_shift;

    // changing peaks height
    change = fRand(0,1) * scen->height_severity;
    new_value = change + peaks_height[j];
    if(new_value < scen->min_height)
      peaks_height[j] = 2.0 * scen->min_height - peaks_height[j] - change;
    else if(new_value > scen->max_height)
      peaks_height[j] = 2.0 * scen->max_height - peaks_height[j] - change;
    else
      peaks_height[j] = new_value;

    // changing peaks width
    change = fRand(0,1) * scen->width_severity;
    new_value = change + peaks_width[j];
    if(new_value < scen->min_width)
      peaks_width[j] = 2.0 * scen->min_width - peaks_width[j] - change;
    else if(new_value > scen->max_width)
      peaks_width[j] = 2.0 * scen->max_width - peaks_width[j] - change;
    else
      peaks_width[j] = new_value;
  }
}

double MovingPeaks::maxValue(std::vector<double> values){
  double result = 0.0;
  for(int i=0; i<scen->npeaks; i++){
    if(values[i]>result)
      result = values[i];
  }
  return result;
}

double MovingPeaks::fRand(double fMin, double fMax){
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

void MovingPeaks::resetProblem() {
  // showPeakes();
  peaks_height.clear();
  peaks_width.clear();
  peaks_position.clear();
  initializePeaks();
}