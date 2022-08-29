/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "GenAlgs.h"
#include <math.h>

using namespace std;

int main (int argc, char *argv[]){

//dichiarazione delle variabili

	int numero_individui = 400;
	int n_generazioni = 1000;
	double elitarieta = 1.5;
	double mutabilita = 0.15;

	int dimensioni = 50;
	int n_riproduzioni = int(numero_individui/4);
	Random rnd;
	Percorso map(&rnd, dimensioni);
	vector<Posizione> pos;
	Posizione citta;
	string state, capital;
	
//inizializzazione delle variabili
	ofstream out ("StatCap.dat");
	fstream in;
  in.open("./American_capitals.dat", ios::in);
  
  if (!in.good()) {
      cerr << "Errore lettura file capitali" << endl;
      return 1;
   }

	for(int i=0; i<50; i++) {
		in >> state >> capital >> citta.x >> citta.y;
		cout << "x: " << citta.x << "	y: " << citta.y << endl;
		pos.push_back(citta);	
	}
	
	in.close();
	
	cout << "numero città: " << pos.size() << endl;
	
	map.SetVec(pos);
	
	Popolazione pop(&rnd, &map, numero_individui, dimensioni, elitarieta, mutabilita);
	
	pop.GeneraPop();
	pop.Statistica();
	cout << pop.GetBest() << "\t\t" << pop.GetMedia()  << "\t\t" << pop.GetStdev() << endl;
	
	for (int i=0; i<n_generazioni; i++){
		pop.NewGen(n_riproduzioni);
		pop.Statistica();

		out << pop.GetBest() << "\t\t" << pop.GetMedia()  << "\t\t" << pop.GetStdev() << endl;
		cout << pop.GetBest() << "\t\t" << pop.GetMedia()  << "\t\t" << pop.GetStdev() << endl;
	}
	
	out.close();
	
	cout << "lunghezza minima: " << pop.GetBest() << endl;	
	map.Stampa("map_cap.dat");
	pop.GetInd(0).Stampa("path_cap.dat");
	
		
	return 0;
}
