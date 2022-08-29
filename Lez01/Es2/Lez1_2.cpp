#ifndef __RANDOM__H

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <cmath>

#include "random.h"
#include "function.cpp"

using namespace std;
 
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
   
   ofstream out1("Unif");
   ofstream out2("Exp");
   ofstream out3("Lor");
   ofstream out;
   
//ESERCIZIO 2: isotgrammare 10^4 somme parziali di ogni distribuzione 
     
//dichiaro le variabili
	//vector<int> N=100;
	vector<double> SnU, SnE, SnL; //alla fine ciascun sn deve contenere 4 elementi
	//sono rispettivamente le somme parziali per le distribuzioni uniforme, esponenziale e Lorentziana
	double M = 10000;
	double N[] = {1,2,10,100};
	int size=4;	
	int lambda=1;
	int Gamma=1;
	int mu=0;
	
	double u, e, l;
	
//Uniforme (standard dice)
	  out.open("Unif");
		for (int i = 0; i < M; i++){
        SnU.resize(0);
        out << i+1;
        //Calcolo 4  versioni diverse per S_N, N = 1,2,10,100
        for(int j = 0; j<size; j++){
            u = 0;
            for(int k = 0; k<N[j]; k++){
                u += rnd.Rannyu();
            }
            SnU.push_back(u/(N[j]) );
            out << "," << SnU.at(j);
        }
        out << endl;
    }
		out.close();
		
		
//Distribuzione esponenziale	
		out.open("Exp");
		for (int i = 0; i < M; i++){
        SnE.resize(0);
        out << i+1;
        //Calcolo 4  versioni diverse per S_N, N = 1,2,10,100
        for(int j = 0; j<size; j++){
            e = 0;
            for(int k = 0; k<N[j]; k++){
                e += rnd.Exp(lambda);
            }
            SnE.push_back(e/(N[j]) );
            out << "," << SnE.at(j);
        }
        out << endl;
    }
		out.close();
		
//Distribuzione Lorentziana
		out.open("Lor");
		for (int i = 0; i < M; i++){
        SnL.resize(0);
        out << i+1;
        //Calcolo 4  versioni diverse per S_N, N = 1,2,10,100
        for(int j = 0; j<size; j++){
            l = 0;
            for(int k = 0; k<N[j]; k++){
                l += rnd.Lorentz(mu, Gamma);
            }
            SnL.push_back(l/(N[j]) );
            out << "," << SnL.at(j);
        }
        out << endl;
    }
		out.close();

  rnd.SaveSeed();
  return 0;
}

#endif
