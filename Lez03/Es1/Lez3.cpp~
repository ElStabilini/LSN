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
#include <algorithm>

#include "random.h"
#include "function.cpp"
#include "Funzioni.h"

using namespace std;

/*Da rivedere: non funziona*/

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
     
/*===========================================================================================
	ESERCIZIO 3 PT.1
	Compute the price of a call or put option price by a direct sampling
	of the price at the expiry.
	
	==========================================================================================
	
	Procedimento
	- sample N times a geometric Brownian motion s(T) at the expiry T;
	- from the sample, estimate the option price as the profit max(S(T)-K,0), K being the 			strike price;
		
=============================================================================================*/
   
//definizione variabili
	double t = 1;
	double S0 = 100; //starting asset price
	double T = 1;	//expiry
	double K = 100; //strike price
	double r = 0.1; //risk free interes rate 
	double sigma = 0.25; //volatility
	double sigma2 = sigma*sigma;
	
	long int M=100000; //numero totale dei lanci              
  long int N=100; //numero di blocchi             
	vector<double> St; //vettore dei numeri campionati secondo la distribuzione finale
	vector<double> CallOpt;
	vector<double> PutOpt; 
	vector<double> x; 
	
	vector<double> aveCall; //medie di ciascun blocco
	vector<double> ave2Call;
	vector<double> sum_progCall; //medie progressive
	vector<double> sum2_progCall;
	vector<double> err_progCall; //errore sui risultati progressivi
	
	vector<double> avePut; //medie di ciascun blocco
	vector<double> ave2Put;
	vector<double> sum_progPut; //medie progressive
	vector<double> sum2_progPut;
	vector<double> err_progPut; //errore sui risultati progressivi

		
	//inizializzo i vettori
	for(int i=0; i<M; i++) {
		St.push_back(rnd.S1(t, S0, r, sigma2)); //riempio r con M numeri casuali
		CallOpt.push_back( exp(-r*T)*max(St.at(i)-K, 0.) );
		PutOpt.push_back( exp(-r*T)*max(K-St.at(i), 0.) );
	}
	
  BlockStat(N, CallOpt, aveCall, ave2Call, sum_progCall, sum2_progCall, err_progCall);
  BlockStat(N, PutOpt, avePut, ave2Put, sum_progPut, sum2_progPut, err_progPut);
  
	Print(sum_progCall, err_progCall, N, "CallOptionPrice_1.dat");
	Print(sum_progPut, err_progPut, N, "PutOptionPrice_1.dat");
	
/*===========================================================================================
	ESERCIZIO 3 PT.2
	Compute the price of a call or put option price by a direct sampling
	of the price at the expiry.
	
	==========================================================================================
	
	Procedimento
	- sample N sequences of a geometric Brownian motion s(t_i) from t_0 = 0 to t_n = N;
	- from the sample, estimate the option price as the profit max(S(T)-K,0), K being the 			strike price;
	
=============================================================================================*/

	aveCall.clear(); //medie di ciascun blocco
	ave2Call.clear();
	sum_progCall.clear(); //medie progressive
	sum2_progCall.clear();
	err_progCall.clear(); //errore sui risultati progressivi
	
	avePut.clear(); //medie di ciascun blocco
	ave2Put.clear();
	sum_progPut.clear(); //medie progressive
	sum2_progPut.clear();
	err_progPut.clear();
	
	int steps = 100;
	double time = T/steps;
		
	for(int i=0; i<M; i++) {
		for(int i=0; i<steps; i++)
			St.push_back(rnd.S2(time, S0, r, sigma2)); //riempio r con M numeri casuali
		CallOpt.push_back( exp(-r*T)*max(St.at(i)-K, 0.) );
		PutOpt.push_back( exp(-r*T)*max(K-St.at(i), 0.) );
	}
	
  BlockStat(N, CallOpt, aveCall, ave2Call, sum_progCall, sum2_progCall, err_progCall);
  BlockStat(N, PutOpt, avePut, ave2Put, sum_progPut, sum2_progPut, err_progPut);
  
//stampo i risultati in colonna separati da una ,
	Print(sum_progCall, err_progCall, N, "CallOptionPrice_2.dat");
	Print(sum_progPut, err_progPut, N, "PutOptionPrice_2.dat");
	
	rnd.SaveSeed();
   
  return 0;
}

#endif
#endif
#endif
