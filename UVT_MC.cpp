#include "UVT_MC.h"

int main()
{
  sp = 10;
  T = 298;
  L = 100;
  int steps = 3;

  create_matter();

  data.open("data.xyz");
  gib_data_bls();

  for (int t = 0; t < steps; t++){
  next_step();
  }
  data.close();
  return 0;
}
