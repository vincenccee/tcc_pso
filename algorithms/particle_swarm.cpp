#include "particle_swarm.hpp"

using namespace std;

ParticleSwarm::ParticleSwarm(Problem *problem, int tamPopulation){
  this->problem = problem;
  this->tamPopulation = tamPopulation;
  this->velocityInertia = 0.6;
  this->acelerationCons1 = 0.20;
  this->acelerationCons1 = 0.20;
}

ParticleSwarm::ParticleSwarm(){}

ParticleSwarm::~ParticleSwarm(){}

void ParticleSwarm::showPopulation(){
  std::vector<double> position;
  for(int i=0; i<tamPopulation; i++){
    cout << "ind " << i << "\t";
    position.clear();
    position = swarm->getIndividual(i)->getCurrentPosition();
    for(unsigned int j=0; j<position.size() ; j++){
      cout << " * " << position[j];
    }
    cout << " - fitness: " << swarm->getIndividual(i)->getFitness() << endl;
  }
}

void ParticleSwarm::evolutionaryCicle(int iterations, int runs){
  this->iterations = iterations;
  this->runs = runs;
  popdata.open("testdata1.txt");
  double result;

  for(int i=0; i<this->iterations; i++){
    bestPopulationFitness.push_back(0.0);
    bestIndividualFitness.push_back(0.0);
    populationDiversity.push_back(0.0);
  }

  for(int j=0; j<this->runs; j++){
    this->swarm = new Population(tamPopulation, problem->getDimension(), problem->getLowerBound(j), problem->getUpperBound(j));
    swarm->initializePopulation();
    this->m_nmdf = 0;
    evaluatePopulationFitness();
    initializeBest();
    cout << "************** inicio ****************" << endl;
    showPopulation();
    for(int i=0; i<this->iterations; i++){
      updateBest(i);
      updatePlot(i);
      updateParticleVelocity();
      updateParticlePosition();
      evaluatePopulationFitness();
    }
    result = this->bestFitness;
    finalFitness.push_back(result);

    cout << "\nBest Finess = " << result << endl;
    problem->resetProblem();
  }
  for(int j=0; j<this->iterations; j++){
    bestPopulationFitness[j] = bestPopulationFitness[j]/this->runs;
    bestIndividualFitness[j] = bestIndividualFitness[j]/this->runs;
    populationDiversity[j] = populationDiversity[j]/this->runs;
  }
  cout << "************** Final ****************" << endl;
  cout << "Media do Fitness: " << arithmeticAverage(finalFitness) << endl;
  cout << "Desvio Padrao: " << standardDeviation(finalFitness) << endl;

  gnu_plot_convergence_best_mean(bestIndividualFitness, bestPopulationFitness, iterations, "MelhoraFitness", "melhora_fit");
  gnu_plot_convergence(populationDiversity, iterations, "pop_diversity", "fitnessdaPopulacao", "Divesidade Genotipica", 1);
  popdata.close();
}

void ParticleSwarm::evaluatePopulationFitness(){
  int newFitness;
  Individual *tmpInd;
  // For para ser paralizado (OPEN-MP)
  // #pragma omp parallel for
  for(int i=0; i<tamPopulation; i++){
    tmpInd = swarm->getIndividual(i);
    newFitness = problem->evaluateFitness(tmpInd->getCurrentPosition());
    if(problem->fitnesIsBetter(newFitness, tmpInd->getFitness())){
      tmpInd->setBestFitness(newFitness);
      tmpInd->setBestPosition(tmpInd->getCurrentPosition());
    }
    tmpInd->setFitness(newFitness);
    swarm->updateIndividual(*tmpInd, i);
  }
}

void ParticleSwarm::updateParticleVelocity(){
  std::vector<double> newVelocity;
  Individual *tmpInd;

  for(int i=0; i<tamPopulation; i++){
    tmpInd = swarm->getIndividual(i);
    newVelocity.clear();
    for(unsigned int j=0; j<problem->getDimension(); j++){
      cout << "bp: " << acelerationCons2 * fRand(0,1) * (this->bestPosition[j] - tmpInd->getPosition(j)) << endl;
      newVelocity.push_back((velocityInertia * tmpInd->getVelocity(j)) +
                            (acelerationCons1 * fRand(0,1) * (tmpInd->getBestPosition(j) - tmpInd->getPosition(j))) +
                            (acelerationCons2 * fRand(0,1) * (this->bestPosition[j] - tmpInd->getPosition(j))));
    }
    swarm->getIndividual(i)->setVelocity(newVelocity);
  }
}

void ParticleSwarm::updateParticlePosition(){
  std::vector<double> newPosition;
  std::vector<double> currentPosition;
  std::vector<double> velocity;
  Individual *tmpInd;
  for(int i=0; i<tamPopulation; i++){
    tmpInd = swarm->getIndividual(i);
    currentPosition = tmpInd->getCurrentPosition();
    newPosition.clear();
    for(int j=0; j<problem->getDimension(); j++){
      newPosition.push_back(currentPosition[j] + tmpInd->getVelocity(j));
    }
    swarm->getIndividual(i)->setCurrentPosition(validatePosition(newPosition));
  }
}

std::vector<double> ParticleSwarm::validatePosition(std::vector<double> position){
  for (unsigned int i = 0; i < position.size(); i++) {
    if(position[i] > problem->getUpperBound(i)){
      position[i] = problem->getUpperBound(i);
    }else if(position[i] < problem->getLowerBound(i)){
      position[i] = problem->getLowerBound(i);
    }
  }
  return position;
}

void ParticleSwarm::initializeBest(){
  this->bestFitness = swarm->getIndividual(0)->getFitness();
  this->bestPosition = swarm->getIndividual(0)->getCurrentPosition();
}

void ParticleSwarm::updateBest(int pos){
  Individual *tmpInd;
  for(int i=0; i < tamPopulation; i++) {
    tmpInd = swarm->getIndividual(i);
    if(problem->fitnesIsBetter(tmpInd->getFitness(), this->bestFitness)){
      this->bestPosition = tmpInd->getCurrentPosition();
      this->bestFitness = tmpInd->getFitness();
    }
  }
  bestIndividualFitness[pos] += this->bestFitness;
}

double ParticleSwarm::fRand(double fMin, double fMax){
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

void ParticleSwarm::updatePlot(int pos){
  double totalFit = 0;
  double mediaFit = 0;

  for(int i=0; i < tamPopulation; i++) {
    totalFit += swarm->getIndividual(i)->getFitness();
  }
  mediaFit = totalFit/tamPopulation;
  bestPopulationFitness[pos] += mediaFit;
  populationDiversity[pos] += defaultGenotypicDiversityMeasure();
}

void ParticleSwarm::gnu_plot_convergence(std::vector<double> mean_gen, int m_gen, std::string name, std::string title, std::string y_axis, double max_range){

  FILE *pipe = popen("gnuplot", "w");

  name = space2underscore(name);
  title = space2underscore(title);

  if (pipe != NULL) {
    fprintf(pipe,"set terminal postscript eps enhanced color colortext font 'Arial,22'\n");
    fprintf(pipe,"set key font 'Arial,18'\n");
    fprintf(pipe,"set encoding iso_8859_1\n");
    fprintf(pipe,"set xlabel 'Gera{/E \347}{/E \365}es'\n");
    string output_y_label("set ylabel '" + y_axis + "'\n");
    fprintf(pipe, "%s", output_y_label.c_str());
    string output("set output '" + OUTPUT_DIR + name + std::string(".eps") + "'\n");
    fprintf(pipe, "%s", output.c_str());

    if(max_range){
      //fprintf(pipe, "set xrange [-1:42]\n");
      string output_range("set yrange [0:" + to_string(max_range) + "]\n");
      fprintf(pipe, "%s", output_range.c_str());
    }
    fprintf(pipe, "set style line 1 lc rgb 'black' \n");

    ofstream color1(OUTPUT_DIR + name + ".dataset");
    for(int k=0; k<m_gen; k++) {
      color1 << k << " " << mean_gen[k] << endl;
    }

    color1.close();
    std::string str_plot = "plot '" + OUTPUT_DIR + name + ".dataset' with lines ls 1 title '" + title + "'\n";
    fprintf(pipe, "%s", str_plot.c_str());

    fflush(pipe);
    pclose(pipe);
  } else{
    std::cout << "Could not open pipe" << std::endl;
  }
}

void ParticleSwarm::gnu_plot_convergence_best_mean(std::vector<double> d_data1, std::vector<double> d_data2, int n_lines, std::string title, std::string filename){

  FILE *pipe = popen("gnuplot", "w");

  title = space2underscore(title);
  filename = space2underscore(filename);

  ofstream melhor_output(OUTPUT_DIR + filename + "_melhor.output");
  ofstream media_output(OUTPUT_DIR + filename + "_media.output");
  for(int k=0; k<n_lines; k++) {
    melhor_output << k << " " << d_data1[k] << endl;
    media_output << k << " " << d_data2[k] << endl;
  }
  melhor_output.close();
  media_output.close();

  if (pipe != NULL) {
    fprintf(pipe, "set bmargin 7\n");
    fprintf(pipe, "unset colorbox\n");
    fprintf(pipe, "set terminal postscript eps enhanced color colortext font 'Arial,22'\n");
    fprintf(pipe, "set key font 'Arial,18'\n");
    fprintf(pipe, "set encoding iso_8859_1\n");
    fprintf(pipe, "set xlabel 'Gera{/E \347}{/E \365}es'\n");
    fprintf(pipe, "set ylabel 'M{/E \351}dia Fitness'\n");
    string output("set output '" + OUTPUT_DIR + filename + std::string(".eps") + "'\n");
    fprintf(pipe, "%s", output.c_str());
    //fprintf(pipe, "set xrange [-1:42]\n"); // set the terminal
    //fprintf(pipe, "set yrange [-1:42]\n"); // set the terminal
    fprintf(pipe, "set style line 1  lc rgb '#B22C2C' dashtype 2 \n");
    fprintf(pipe, "set style line 2  lc rgb '#0404E9' lt 2 \n");
    std::string str_title = "set title '" + title + "'\n";
    fprintf(pipe, "%s", str_title.c_str());
    std::string str_plot1 = "plot '" + OUTPUT_DIR + filename + "_melhor.output' using 1:2 title 'Melhor' ls 1 lw 5 with lines, ";
    std::string str_plot2 = "'" + OUTPUT_DIR + filename + "_media.output' using 1:2 title 'M{/E \351}dia' ls 2  lw 3 with lines\n ";
    fprintf(pipe, "%s", str_plot1.c_str());
    fprintf(pipe, "%s", str_plot2.c_str());

    fflush(pipe);
    pclose(pipe);
  } else{
    std::cout << "Could not open pipe" << std::endl;
  }
}

std::string ParticleSwarm::space2underscore(std::string text) {
  for(std::string::iterator it = text.begin(); it != text.end(); ++it) {
      if(*it == ' ') {
          *it = '_';
      }
  }
  return text;
}

double ParticleSwarm::defaultGenotypicDiversityMeasure(){
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

double ParticleSwarm::standardDeviation(std::vector<double> data){
  double media = arithmeticAverage(data);
  double sum1 = 0;

  for (unsigned int i=0; i<data.size(); i++) {
      sum1 += pow((data[i] - media), 2);
  }
  return sqrt(sum1 / (double)(data.size() - 1));
}

double ParticleSwarm::arithmeticAverage(std::vector<double> data){
  double sum = 0;
  for (unsigned int i=0; i<data.size(); i++) {
      sum += data[i];
  }
  return (sum / data.size());
}