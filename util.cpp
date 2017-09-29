#include "util.hpp"

using namespace std;

Util::Util(){}
Util::~Util(){}

void Util::openData(string filename){
  popdata.open(filename);
}
void Util::closeData(){
  popdata.close();
}

double Util::fRand(double fMin, double fMax){
  double number;
  if(RANDOM_DISTRIBUTION == 1) {
    double f = (double)rand() / RAND_MAX;
    number = fMin + f * (fMax - fMin);
  } else {
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_real_distribution<double> uni(fMin,fMax); // guaranteed unbiased
    number = uni(rng);
  }
  return number;
}

void Util::gnu_plot_convergence(std::vector<double> mean_gen, int m_gen, std::string name, std::string title, std::string y_axis, double max_range){

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

void Util::gnu_plot_convergence_best_mean(std::vector<double> d_data1, std::vector<double> d_data2, int n_lines, std::string title, std::string filename){

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

std::string Util::space2underscore(std::string text) {
  for(std::string::iterator it = text.begin(); it != text.end(); ++it) {
      if(*it == ' ') {
          *it = '_';
      }
  }
  return text;
}

double Util::standardDeviation(std::vector<double> data){
  double media = arithmeticAverage(data);
  double sum1 = 0;

  for (unsigned int i=0; i<data.size(); i++) {
      sum1 += pow((data[i] - media), 2);
  }
  return sqrt(sum1 / (double)(data.size() - 1));
}

double Util::arithmeticAverage(std::vector<double> data){
  double sum = 0;
  for (unsigned int i=0; i<data.size(); i++) {
      sum += data[i];
  }
  return (sum / data.size());
}