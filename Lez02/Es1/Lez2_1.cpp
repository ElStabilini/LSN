#ifndef __RANDOM__H
#define __RANDOM__H

#ifndef __INTEGRAZIONE_H__
#define __INTEGRAZIONE_H__

#ifndef __FUNZIONI__H
#define __FUNZIONI__H

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <cmath>

#include "random.h"
#include "Funzioni.h"
#include "Integrazione.h" 
#include "function.cpp"

using namespace std;

bool Lancio(double x1, double x2, double y2);
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
 
//DEFINIZIONE DELLE VARIABILI
	long int M=500000; //numero totale di stime dell'integrale            
  long int N=500; //numero di blocchi             
	vector<double> x; //vettore indici
	
	vector<double> media; //vettore in cui raccolgo tutte le stime di I con metodo media
	vector<double> aveM; //vettore di N valori - contiene le medie di ciascun blocco
	vector<double> ave2M; //vettore con gli stessi elementi ma al quadrato
	vector<double> sum_progM; //vettore di N valori - contiene le medie progressive (ciascun elemento progressiva fino ad N-1)
	vector<double> sum2_progM; //vettore con gli stessi elementi ma al quadrato
	vector<double> err_progM; //errore sui risultati progressivi
	
	vector<double> ImpSamp; //raccolgo tutte le stime di I con importance sampling
	vector<double> aveI; 
	vector<double> ave2I;
	vector<double> sum_progI; 
	vector<double> sum2_progI; 
	vector<double> err_progI;
	
	
	for(int i=0; i<N; i++) 
		x.push_back(i); //creo il vettore [0,...,N-1]
   
//PRIMA PARTE: Calcolo dell'integrale con il metodo della media con 1000 punti

	FunzioneBase* f = new Coseno(M_PI/2, M_PI/2, 0, 0);
	Integral myintegral;
	
	for(int i=0; i<M; i++)
			media.push_back(myintegral.Media(f, 0, 1, 1000, rnd));
			
	BlochStat(N, media, aveM, ave2M, sum_progM, sum2_progM, err_progM);
	Print(sum_progM, err_progM, N, "Media.dat");

//SECONDA PARTE: calcolo dell'integrale con il metodo di importance sampling


	FunzioneBase* pdf = new Retta(-2, 2); 
	
	for(int i=0; i<M; i++)
			ImpSamp.push_back(myintegral.ImpSampling(f, pdf, 0, 1, 1000, rnd));
			
	BlochStat(N, ImpSamp, aveI, ave2I, sum_progI, sum2_progI, err_progI);
	Print(sum_progI, err_progI, N, "ImpSampling.dat");

	rnd.SaveSeed();
	
	//
   
  return 0;
}

#endif
#endif
#endif
