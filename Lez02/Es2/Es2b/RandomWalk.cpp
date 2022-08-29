#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>

#include <vector>
#include <string>

#include "random.h"
#include "position.h"

using namespace std;

//RANDOM WALK CONTINUO con length-step e in dimensione "dim"
vector<double> WalkCon(int length, Random& rnd) {

	double R=1; //lunghezza dello step
	Position pos0(0, 0, 0);
	vector<Position> traj; 
	
//in questo caso, dato che mi muovo nel continuo ho bisogno di estrarre almeno una coppia di numeri per definire la direzione
	
	traj.push_back(pos0);
	
	for(int i=0; i<length; i++) {
		double theta = acos(1-2*rnd.Rannyu(0.,1.));
		double phi = 2*M_PI*rnd.Rannyu(0.,1.);
		double newX = traj.at(i).getX() + R*sin(theta)*cos(phi);
		double newY = traj.at(i).getY() + R*sin(theta)*sin(phi);
		double newZ = traj.at(i).getZ() + R*cos(theta);

		traj.push_back( Position(newX, newY, newZ) );
	}

	vector<double> distance;
	for(int i=0; i<length; i++)
		distance.push_back( traj.at(i).Distanza(pos0) );
	
	return distance;
};


/*
vector<double> x (length, rnd.Rannyu(-R,R));
	vector<double> y;
	vector<double> z;
	
	traj.push_back(pos0);
	
	for(int i=0; i<length; i++){
		double ri = sqrt(R*R - pow(x.at(i),2));
		y.push_back(rnd.Rannyu(-ri,ri));
		z.push_back(rnd.Rannyu()*sqrt(ri*ri - pow(y.at(i),2)));

		double newX = traj.at(i).getX() + x.at(i);
		double newY = traj.at(i).getY() + y.at(i);
		double newZ = traj.at(i).getZ() + z.at(i);
		traj.push_back( Position(newX, newY, newZ) );
	}*/


