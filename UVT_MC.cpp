#include "UVT_MC.h"
using namespace std;

int main(int argc, char **argv)
{
  if (argc != 2){
    cout << "Please include a file of conditions." << endl;
    return 1;
  }
  srand(time(NULL));
  system_t system;
  vector<double> conditions;
  double temp;
  ifstream input (argv[1]);
  while (input >> temp){
    conditions.push_back(temp);
  }
  system.sp = conditions[0];
  system.L = conditions[1];
  system.T = conditions[2];
  system.steps = conditions [3];

  ofstream data;
  create_matter(&system);
  data.open("data.xyz");
  gib_data_bls(&system, data);
  for (int t = 0; t < system.steps; t++){
    next_step(&system, data);
  }
  data.close();
  return 0;
}
