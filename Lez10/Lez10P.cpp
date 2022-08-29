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
#include "randomP.h"
#include "GenAlgs.h"
#include <mpi.h>
#include <math.h>

#define ROOT1 0
#define ROOT2 1
#define ROOT3 2
#define ROOT4 3

using namespace std;

int main (int argc, char *argv[]){

	int size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//dichiarazione delle variabili
	int numero_individui = 400;
	int n_generazioni = 1000;
	double select = 1.5;
	double mut = 0.15;

	int scambio = 40; //questo è il numero di generazioni dopo cui avviene lo scambio
	
	//"migrazioni" che avvengono tra i continenti, ovvero il numero di volte per cui i diversi nodi si scambiano gli individui migliori della popolazione
	int dimensioni = 50;
	int n_riproduzioni = int(numero_individui/4);
	Random rnd(rank);
	Percorso map(&rnd, dimensioni);
	vector<Posizione> pos;
	Posizione citta;
	string state, capital;
	
	ofstream out ("CapBest"+to_string(rank)+".dat");
	
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
	map.Stampa("map_capP.dat");
	Popolazione pop(&rnd, &map, numero_individui, dimensioni, select, mut);
	
	//Da qui in avanti il codice è parallelizzato perchè in ogni continente (su ogni nodo) la popolazione evolve in maniera indipendente
	
	pop.GeneraPop();
	pop.Statistica();
	cout << pop.GetBest() << "\t\t" << pop.GetMedia()  << "\t\t" << pop.GetStdev() << endl;
	
	//in questo caso in pratica conto 25 volte fino a 40
	for (int i=0; i<n_generazioni/scambio; i++){
		for(int j=0; j<scambio; j++) { //per 40 generazioni evolvono indipendentemente
			pop.NewGen(n_riproduzioni);
			pop.Statistica();

			out << pop.GetBest() << "\t\t" << pop.GetMedia()  << "\t\t" << pop.GetStdev() << endl;
		}
		
		//creo quattro matrici una per ogni rank, ciascuna matrice numero_individui righe e dimensioni colonne (accosto in verticale i DNA dei diversi individui
		int zer[numero_individui][dimensioni];
		int uno[numero_individui][dimensioni];
		int due[numero_individui][dimensioni];
		int tre[numero_individui][dimensioni];

		//inizializzo tutte le matrici
		for (int k=0; k<numero_individui; k++) {
			for(int l=0; l<dimensioni; l++) {
				zer[k][l] = pop.GetInd(k).GetGene(l);
        uno[k][l] = pop.GetInd(k).GetGene(l);
        due[k][l] = pop.GetInd(k).GetGene(l);
        tre[k][l] = pop.GetInd(k).GetGene(l);
			}
		}
		
		MPI_Bcast(zer, dimensioni*numero_individui, MPI_INT, ROOT1 , MPI_COMM_WORLD);
		MPI_Bcast(uno, dimensioni*numero_individui, MPI_INT, ROOT2 , MPI_COMM_WORLD);
		MPI_Bcast(due, dimensioni*numero_individui, MPI_INT, ROOT3, MPI_COMM_WORLD);
		MPI_Bcast(tre, dimensioni*numero_individui, MPI_INT, ROOT4, MPI_COMM_WORLD);
		
//syntax: MPI_BCAST(elementi da trasferire, numero di parametri da trasferire, datatype, root, communicator)
//syntax: MPI_BCAST(buffer, int count, MPI_Datatype, root, communicator)

		int tutti[3][numero_individui][dimensioni];
//accosto le matrici una all'altra (creo una specie di tensore)
		for (int k=0; k<numero_individui; k++) {
			for(int l=0; l<dimensioni; l++) {
      	tutti[0][k][l] = zer[k][l];
				tutti[1][k][l] = uno[k][l];
        tutti[2][k][l] = due[k][l];
        tutti[3][k][l] = tre[k][l];       
      }
		}
		// ognuno riarrangia e prende i suoi nuovi

		for (int k=0; k<numero_individui; k++){
			vector<int> vec; 
			vec.empty(); //svuoto dal ciclo precedente
			int a = (k+rank)%4; //mi assicuro che nessuno ripeschi se stessp
				for(int l=0; l<dimensioni; l++) //creo il DNA del nuovo individuo
					vec.push_back(tutti[a][k][l]);
			
			Individuo nuovo(&rnd,dimensioni,vec); //creo un nuovo individuo
			nuovo.Fit(pop.GetMap());
			pop.Sostituisci(k, nuovo);
		}
	}
	

	
	cout << "lunghezza minima percorso: " << pop.GetBest() << endl;	
	pop.GetInd(0).Stampa("path_capP.dat");
	
	MPI_Finalize();
	out.close();	
		
	return 0;
}
