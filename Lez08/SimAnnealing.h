#define PREC 0.0015
#define delta 2

using namespace std;

ofstream out; 

double accepted_ann, accepted_met, attempted_ann, attempted_met;
double muold, sigmaold;
double m=1, hbar=1;
double beta, delta_ann;
int eq_step;
double mu = 0.5;
double sigma = 0.5;


//funzioni
void InitializeRND(Random& rnd);
double Psi2(double x);
double Psi(double x);
double V(double x);
double derivative(double x);
double K(double x);
double Eloc(double x);

double Metropolis(Random& rnd, double x);
double Estrazione(Random& rnd, double& x);
double Integral(Random& rnd, int M, int L, bool print);
double Evol(Random &rnd, double beta);
