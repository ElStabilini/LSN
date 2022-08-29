#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>

#include <vector>
#include <string>

#include "random.h"
#include "function.cpp"
#include "RandomWalk.cpp"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;

   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
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


	//definizione variabili   
  long int M = 100000; //numero di RW             
  long int N = 100; //numero di blocchi             
  int length = 100;
  
  vector<double> ave;
  vector<double> ave2;
  vector<double> sum_prog;
  vector<double> sum2_prog;
  vector<double> err_prog;
  vector< vector<double> > RW; //vettore di M elementi
  vector< vector<double> > data; //vettore di length elementi
  vector<double> distanza;
  vector<double> err;


//genero tutti i RW   
  for(int i=0; i<M; i++)
  	RW.push_back(WalkCon(length, rnd));
  
//ordino la matrice raggruppanndo per step
  for(int l=0; l<length; l++){
  
  	vector<double> step;
  	for(int j=0; j<M; j++)
  		step.push_back(RW.at(j).at(l));
  		
  	data.push_back(step); //prendo lo stesso step (l) per tutti i RW (j)
  }
 
 	for(int i=0; i<length; i++) {
  	BlockStat(N, data.at(i), ave, ave2, sum_prog, sum2_prog, err_prog);
		distanza.push_back(sum_prog.at(i));
		err.push_back(err_prog.at(i));
	}
	
	

	Print(distanza, err, length, "RWR3.dat");

  rnd.SaveSeed();
  return 0;
}
