#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "SimAnnealing.h"

#define _USE_MATH_DEFINES

using namespace std;

int main (int argc, char *argv[]){

  Random rnd;
	InitializeRND(rnd);

  double H;
  
  int L = 15; // Numero di punti estratti in ogni blocco (a temperatura costante)
  vector <double> mus, sigmas, betas; 
  double betamin = 1;
  double dbeta = 2.5;
  double err=0.5;
  int k=0, iter=0;

  out.open("annealing.dat");
  double error = PREC,  beta;
  double aveH=0, ave2H=0;


  while(err >= PREC) {
  	iter++;
    
    betas.push_back(betamin + k*dbeta); // Aggiorno la temperatura
    beta = betas.at(k);
    delta_ann = 1./sqrt(beta);
    double sum_H=0, sum2_H=0, sum_mu=0, sum_sigma=0;
    
 
  //Ciclo nel blocco a temperatura costante
    for(int j = 0; j < L; j++){
      H = Evol(rnd, beta);
      sum_H += H;
      sum2_H += H*H;
      sum_mu += mu;
      sum_sigma += sigma;
    }
    
   	// Deviazione standard sulla stima di H  a temperatura costante.
    err = sqrt(sum2_H/L - sum_H*sum_H/L/L);
    sum_H = sum_H/L;
    
    //stime progressive
    aveH += sum_H;
    ave2H += sum_H*sum_H;
    error = k ? sqrt((ave2H - pow(ave2H,2))/k) : PREC;
    
    // Stime mu e sigma
    mus.push_back(sum_mu/L);
    sigmas.push_back(sum_sigma/L);

  // Formato: #blocco, media blocco, errore blocco
    out << k+1 << "," << aveH/(k+1) << "," << error << "," << sum_H << "," << err << "," << mus.at(k) << "," << sigmas.at(k) << endl;

    cout << endl << "BLOCCO " << k << ", t = " << 1./beta << ", passo = " << delta_ann << endl;
    cout << " mu = " << mu << ", sigma = " << sigma << endl;
    cout << " Tasso accettazione metropolis: " << (double)accepted_met/attempted_met << endl;
    
    k++;	
		cout << "==============================================================" << endl;
		
		if(iter==1e4) { //inserire controllo su iterazioni
			cout << "unable to find result with precision " << PREC << endl;
			break;
		}
  }

  cout << endl << "mu = " << mu << endl << "sigma = " << sigma << endl;
  out.close();

//PUNTO 3: mostro la variazione della stima di <H> in funzione dei blocchi
	out.open("Punto3.dat");
  double Hfixed = Integral(rnd, 100000, 100, true);
  
  

//PUNTO 4: preparazione dell'istogramma per mostrare la densità di probabilità campionata
  out.open("histo.dat");
  mu = mus.at(mus.size()-1);
  sigma = sigmas.at(sigmas.size() - 1);
  int max = 50000;
  double x = 1;
  for(int i = 0; i < max; i++){
    x = Metropolis(rnd, x);
    out << Metropolis(rnd, x) << endl;
  }

  out.close();
  
  rnd.SaveSeed();
  return 0;
}


/**********************************************************************************************

IMPLEMENTAZIONE FUNZIONI

***********************************************************************************************/


void InitializeRND(Random& rnd){
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
}

double Psi(double x) {
	double f1 = (x-mu)/(M_SQRT2*sigma);
	double f2 = (x+mu)/(M_SQRT2*sigma);
	double f = exp(-f1*f1) + exp(-f2*f2);
	return f;
};

double Psi2(double x) {
	double f1 = (x-mu)/(M_SQRT2*sigma);
	double f2 = (x+mu)/(M_SQRT2*sigma);
	double f = exp(-f1*f1) + exp(-f2*f2);
	return f*f;
};

double V(double x) {
	return pow(x, 4) - 5*x*x/2;
};

double derivative(double x) {
  double alpha = pow((x+mu)/sigma,2);
  double beta = pow((x-mu)/sigma,2);
	return -1./sigma/sigma * ( alpha*exp(-alpha/2) + beta*exp(-beta/2) - exp(-alpha/2) - exp(-beta/2) ); 
};

double K(double x) {
	return derivative(x)*hbar*hbar/(2*m*Psi(x));
};

double Eloc(double x) {
	return V(x) + K(x);
};

double Metropolis(Random& rnd, double x) {
	attempted_met++;
	double x_new = x + rnd.Rannyu(-delta,delta);
  double p = fmin(1., Psi2(x_new)/Psi2(x) ); 
  double y = rnd.Rannyu();
  
  if(y < p){
    accepted_met++;
    x = x_new;
  }
  
  return x;
};

double Estrazione(Random& rnd, double& x){
    x = Metropolis(rnd, x);    
    return Eloc(x);
};


double Integral(Random &rnd, int M, int L, bool print) {
  int N = int(M/L);
  
  eq_step = 500;
  double x = 1, sum; // Punto del campionamento con metropolis
  double inc; //incertezza: errore progressivo

  // Equilibrazione
  for(int i=0; i<eq_step; i++)
    Estrazione(rnd, x);
  
  double avg=0;
  for(int k = 0; k<N; k++) {
  
    sum = 0;
    for(int j = 0; j<L; j++)
      sum += Estrazione(rnd, x);
  
    sum = sum/L;
    avg += sum;

    if(print) {
      double avg2;
      avg2 += sum*sum;
			inc = k ? sqrt((avg2/(k+1)-pow(avg/(k+1),2))/k) : 0;
      out << k+1 << "," << avg/(k+1) << "," << inc << endl;
    } 
  }

  return avg/N;
};


double Evol(Random &rnd, double beta){
  attempted_ann++;
  
	muold = mu;
	sigmaold = sigma;
	double H = Integral(rnd, 70000, 1, false);
	
	mu += rnd.Rannyu(-delta_ann, delta_ann);
	sigma += rnd.Rannyu(-delta_ann, delta_ann);
	double Hnew = Integral(rnd, 70000, 1, false);
  
  double P = (Hnew-H > 0) ? exp(-beta * (Hnew-H)) : 1;
  double y = rnd.Rannyu();

	if(Hnew-H > 0)
  	
  	if(y < P){
    	accepted_ann ++;
    	H = Hnew;
  	}
  	else {
    	mu = muold;
    	sigma = sigmaold;
  	}
  
  return H;
};
