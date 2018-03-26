#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <vector>

struct particle{
  double x[3];
  std::string type;
};

enum stepType {move, add, destroy};

typedef struct{

  double L,
         T,
         k = 0.00138;
  int sp,
      steps;

  double sigma,
         epsilon;

  double pe;

  double dx,
         dy,
         dz;

  bool good;

  stepType step;

  std::vector<particle> particles;
  int p;
  double x,
         y,
         z;
  int choose;
} system_t;

std::ofstream energies;

double get_distance(system_t * sys, int a, int b);
void gib_data_bls(system_t * sys, std::ofstream & data);
double get_pe(system_t * sys);
void evaluate_pe(system_t * sys, double new_pe, double old_pe);
void no_leave_box(system_t * sys);

double get_random_number(int a, int b)
{
  return a + (random()/((double)RAND_MAX)*(b - a));
}

void create_matter(system_t * sys)
{
  for (int i = 0; i < sys->sp; i++){
    particle temp;
    sys->particles.push_back(temp);
  }
  for (int n = 0; n < sys->sp; n++){
    sys->particles[n].x[0] = get_random_number(0, sys->L);
    sys->particles[n].x[1] = get_random_number(0, sys->L);
    sys->particles[n].x[2] = get_random_number(0, sys->L);
  }
}

void move_particle(system_t * sys)
{
  sys->dx = get_random_number(0, 0.5*sys->L),
    sys->dy = get_random_number(0, 0.5*sys->L),
    sys->dz = get_random_number(0, 0.5*sys->L);

  sys->p = rand() % sys->particles.size();

    sys->choose = rand() % 2;

  if (sys->choose == 0){
    sys->particles[sys->p].x[0] += sys->dx;
    sys->particles[sys->p].x[1] += sys->dy;
    sys->particles[sys->p].x[2] += sys->dz;
  }
  else {
    sys->particles[sys->p].x[0] -= sys->dx;
    sys->particles[sys->p].x[1] -= sys->dy;
    sys->particles[sys->p].x[2] -= sys->dz;
  }
  no_leave_box(sys);
}

void no_leave_box(system_t * sys)
{
  int np = sys->particles.size();

  for (int n = 0; n < np; n++){
    for (int i = 0; i <= 2; i++){
      if (sys->particles[n].x[i] > sys->L){
        sys->particles[n].x[i] -= sys->L;
      }
      if (sys->particles[n].x[i] < 0){
        sys->particles[n].x[i] += sys->L;
      }
    }
  }
}

void add_particle(system_t* sys)
{
  particle temp;
  sys->particles.push_back(temp);

  sys->particles.back().x[0] = get_random_number(0, sys->L);
  sys->particles.back().x[1] = get_random_number(0, sys->L);
  sys->particles.back().x[2] = get_random_number(0, sys->L);
}

void remove_particle(system_t * sys)
{
  if(!sys->particles.size()){ return;}
  sys->p = rand() % sys->particles.size();
  sys->x = sys->particles[sys->p].x[0];
  sys->y = sys->particles[sys->p].x[1];
  sys->z = sys->particles[sys->p].x[2];
  sys->particles.erase(sys->particles.begin() + sys->p);
}

void next_step(system_t * sys, std::ofstream & data)
{
  double old_pe = get_pe(sys);
  sys->particles = sys->particles;

  double pick;
  if(sys->particles.size() == 1) {
    pick = get_random_number(0,2);//avoids floating point errror
  }
  else {
    pick = get_random_number(0, 3);
  }

  if (pick < 1){
    add_particle(sys);
    sys->step = add;
  }
  else if ((1 < pick) && (pick < 2)){
    move_particle(sys);

    sys->step = move;
  }
  else {
    remove_particle(sys);
    sys->step = destroy;
  }

  double new_pe = get_pe(sys);
  evaluate_pe(sys, new_pe, old_pe);

  if (sys->good == false){ 
    if (sys->step == add){
      std::cout << "wowie" << std::endl;
      sys->particles.erase(sys->particles.end()-1);
    }
    else if (sys->step == move){
      std::cout << "hi" << std::endl;
      if (sys->choose == 0){
        sys->particles[sys->p].x[0] -= sys->dx;
        sys->particles[sys->p].x[1] -= sys->dy;
        sys->particles[sys->p].x[2] -= sys->dz;
      }
      else {
        sys->particles[sys->p].x[0] += sys->dx;
        sys->particles[sys->p].x[1] += sys->dy;
        sys->particles[sys->p].x[2] += sys->dz;
      }
    }
    else if (sys->step == destroy){
      std::cout << "hello" << std::endl;
      particle temp;
      temp.x[0] = sys->x;
      temp.x[1] = sys->y;
      temp.x[2] = sys->z;
      sys->particles.push_back(temp);
    }
    else{
      std::cout << "uh oh" << std::endl;
    }
  }
  gib_data_bls(sys, data);
}

double get_distance(system_t * sys, int a, int b, const std::vector<particle> & v)
{
  double d;
  double change_x;
  double change_y;
  double change_z;

  change_x = fabs(v[a].x[0] - v[b].x[0]);
  if (change_x > 0.5*sys->L) {
    change_x -= 0.5*sys->L;
  }
  change_y = fabs(v[a].x[1] - v[b].x[1]);
  if (change_y > 0.5*sys->L) {
    change_y -= 0.5*sys->L;
  }
  change_z = fabs(v[a].x[2] - v[b].x[2]); 
  if (change_z > 0.5*sys->L) {
    change_z -= 0.5*sys->L;
  }
  double change2_x = change_x*change_x;
  double change2_y = change_y*change_y;
  double change2_z = change_z*change_z;

  d = sqrt(change2_x + change2_y + change2_z);
  return d;
}

double get_pe(system_t * sys)
{
  sys->sigma = 0.03345;
  double s = sys->sigma,
         s2 = s*s,
         s6 = s2*s2*s2,
         s12 = s6*s6;

  sys->epsilon = 0.173;
  double e = sys->epsilon;
  sys->pe = 0;

  int np = sys->particles.size();

  for (int a = 0; a < np; a++){
    for (int b = a + 1; b < np; b++){
      double r = get_distance(sys, a, b, sys->particles),
             r2 = r*r,
             r6 = r2*r2*r2,
             r12 = r6*r6;
      sys->pe += 4*e*((s12/r12) - (s6/r6));
    }
  }
  return sys->pe;
}

void evaluate_pe(system_t * sys, double new_pe, double old_pe)
{
  double change_pe = new_pe - old_pe;
  double prob = exp(-1*change_pe*(1/(sys->k*sys->T)));
  if (prob > ((double)random()/(double)RAND_MAX)){
    sys->good = true;
  }
  else {
    sys->good = false;
  }
  std::cout << prob << std::endl;
}

void gib_data_bls(system_t * sys, std::ofstream & data)
{
  int np = sys->particles.size();
  data << np << "\n\n";

  for (int n = 0; n < np; n++){
    data << "Ar " << sys->particles[n].x[0] << " " <<
      sys->particles[n].x[1] << " " << sys->particles[n].x[2] << sys->good 
      << "\n";
  }
}

