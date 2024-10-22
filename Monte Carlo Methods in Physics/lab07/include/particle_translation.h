#ifndef PARTICLE_TRANSLATION_H
#define PARTICLE_TRANSLATION_H

#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<cmath>


using namespace std;


void vector_rotation(double tx, double ty, double tz, double & x, double & y, double & z, double teta);

void particle_translation(double & x1,double & y1, double & x2,double & y2, 
				  double xr,double yr, double Rr, double xa,double ya, double Ra, 
				  int & exist ,double & length);

#endif // PARTICLE_TRANSLATION_H