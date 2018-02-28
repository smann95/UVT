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

double L;
double T;
int sp;
double k = 1;

std::ofstream data;

std::vector<particle> particles;
std::vector<particle> particles_projection;

double get_distance(int a, int b, const std::vector<particle> & v);
void gib_data_bls();
double get_pe();
bool evaluate_pe(double old_pe, double new_pe);

double get_random_number(int a, int b)
{
  return a + (random()/((double)RAND_MAX)*(b - a));
}

void create_matter()
{
  for (int i = 0; i < sp; i++){
    particle temp;
    particles.push_back(temp);
  }
  for (int n = 0; n < sp; n++){
    particles[n].x[0] = get_random_number(0, L);
    particles[n].x[1] = get_random_number(0, L);
    particles[n].x[2] = get_random_number(0, L);
  }
}

void move_particle()
{
  double dx = get_random_number(0, 0.5*L),
         dy = get_random_number(0, 0.5*L),
         dz = get_random_number(0, 0.5*L);

  int n = rand() % particles_projection.size();

  int r = rand() % 2;

  if (r == 0){
    particles_projection[n].x[0] += dx;
    particles_projection[n].x[1] += dy;
    particles_projection[n].x[2] += dz;
  }
  else {
    particles_projection[n].x[0] -= dx;
    particles_projection[n].x[1] -= dy;
    particles_projection[n].x[2] -= dz;
  }
}

void add_particle()
{
  particle temp;
  particles_projection.push_back(temp);

  particles_projection.back().x[0] = get_random_number(0, L);
  particles_projection.back().x[1] = get_random_number(0, L);
  particles_projection.back().x[1] = get_random_number(0, L);
}

void remove_particle()
{
  int n = rand() % particles_projection.size();

  particles_projection.erase(particles_projection.begin() + n);
}

void next_step()
{
 double old_pe = get_pe();

 particles_projection = particles;

  double pick = get_random_number(1, 3);

  if (pick < 2){
    add_particle();
  }
  if ((1 < pick) && (pick < 2)){
    move_particle();
  }
  if (pick <= 3){
    remove_particle();
  }

 double new_pe = get_pe();

  if (evaluate_pe(old_pe, new_pe) == true){
    particles = particles_projection;
  }
  gib_data_bls();
}
double get_distance(int a, int b, const std::vector<particle> & v)
{
  double d;
  double change_x;
  double change_y;
  double change_z;

  change_x = fabs(&v[a].x[0] - &v[b].x[0]);
  if (change_x > 0.5*L) {
    change_x -= 0.5*L;
  }
  change_y = fabs(&v[a].x[1] - &v[b].x[1]);
  if (change_y > 0.5*L) {
    change_y -= 0.5*L;
  }
  change_z = fabs(&v[a].x[2] - &v[b].x[2]); 
  if (change_z > 0.5*L) {
    change_z -= 0.5*L;
  }
  double change2_x = change_x * change_x;
  double change2_y = change_y * change_y;
  double change2_z = change_z * change_z;

  d = sqrt(change2_x + change2_y + change2_z);
  return d;
}

double get_pe()
{
  double sigma = 3.345,
         s = sigma,
         s2 = s*s,
         s6 = s2*s2*s2,
         s12 = s6*s6;

  double epsilon = 1.73e-21,
         e = epsilon;
  double pe = 0;

  int np = particles_projection.size();

  for (int a = 0; a < np; a++){
    for (int b = a + 1; b < np; b++){
      double r = get_distance(a, b, particles_projection),
             r2 = r*r,
             r6 = r2*r2*r2,
             r12 = r6*r6;
      pe += 4*e*((s12/r12) - (s6/r6));
    }
  }
  return pe;
}

bool evaluate_pe(double old_pe, double new_pe)
{
  bool good;
  double change_pe = new_pe - old_pe;
  double prob = exp(-1*change_pe*(1/(k*T)));
  if (prob > ((double)random()/(double)RAND_MAX)){
    good = true;
  }
  else {
    good = false;
  }
  return good;
}

void gib_data_bls()
{
  int np = particles.size();
  data << np << "\n\n";

  for (int n = 0; n < np; n++){
    data << "Ar " << particles[n].x[0] << " " <<
      particles[n].x[1] << " " << particles[n].x[2] << "\n";
  }
}
