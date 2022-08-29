#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>

#include <vector>
#include <string>

#include "random.h"
#include "position.h"

using namespace std;

//RANDOM WALK SU RETICOLO con length-step e in dimensione "dim"
vector<double> WalkRet(int length, Random& rnd) {

	double l=1; //lunghezza dello step
	//Random rnd;
	Position pos0(0, 0, 0);
	vector<Position> traj; 
	vector<int> axis;
	vector<int> sign;
	
	for(int i=0; i<length; i++) {
		axis.push_back(rnd.Int(1,4));
		sign.push_back(rnd.Choice(-1,1));
	}
	
	traj.push_back(pos0);
	
	for(int i=0; i<length; i++) {
	
		//coordinate del posto da cui parte il passo
		double x = traj[i].getX();
		double y = traj[i].getY();
		double z = traj[i].getZ();
		
		/*eg. al passo 0 aggiungo un elemento al vettore (che adesso ha indici 0 e 1, 
			ovvero i e i+1, in generale avrà i+1 elementi)
			questo elemento è una copia del precedente
			evolvo questo elemento (quello con indice i+1) portandomi in una nuova poszione)*/
		
		
		// NB: inserisci un modo per accedere ai metodi privati e modificarli
		if (axis[i] == 1) {
			double newX = x + sign[i]*l;
			traj.push_back( Position(newX, y, z) );
		} else if (axis[i] == 2) {
			double newY = y + sign[i]*l;
			traj.push_back( Position(x, newY, z) );
		} else if (axis[i] == 3) {
			double newZ = z + sign[i]*l;
			traj.push_back( Position(x, y, newZ) );
		}
		else
			cout << "errore, step non intero" << endl; 
	}
	
	vector<double> distance;
	for(int i=0; i<length; i++)
		distance.push_back( traj.at(i).Distanza(pos0) );
	
	return distance;
};
