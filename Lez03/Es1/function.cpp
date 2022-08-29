#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <cmath>

#include "random.h"

using namespace std;

// =========================================================================
//FUNZIONI STATISTICHE
//==========================================================================


//Calcola MEDIA
template <typename T> double CalcolaMedia( const vector<T> & v) {
double sum = 0;
	for (int i=0; i<v.size(); i++)
		sum = sum + v[i];
	return sum/v.size();
};


//Calcola MEDIANA
template <typename T> double CalcolaMediana( vector<T> v ) {
	
	sort( v.begin(), v.end() ) ;  
  double mediana = 0;
	double x = v.size()/2;
	if (v.size()%2 == 0) {
		mediana = (v[v.size()/2]+v[v.size()/2 -1])/2;
	} else {
		int n = (int)x;
		mediana = v[x];
	}
  return mediana;
	
};


//Calcola VARIANZA
template <typename T> double CalcolaVarianza( const vector<T> & v) {
	double scarto = 0;
	double media = CalcolaMedia(v);
	for (int i=0; i<v.size(); i++)
  	scarto = scarto + pow(v[i] - media, 2);
	return scarto/v.size();
};


//Calcola SIGMA
template <typename T> double CalcolaSigma( const vector<T> & v) {
	double mediaq = pow(CalcolaMedia(v), 2);
	vector<T> quadrato;
	for (int i=0; i<v.size(); i++)
		quadrato.push_back(pow(v[i], 2));
	double qmedia = CalcolaMedia(quadrato);
	return sqrt(qmedia - mediaq);
};


//Calcola CORRELAZIONE
template <typename T> double CalcolaCorrelazione (const vector<T> & x, const vector<T> & y) {
	vector<T> xy;
	for (int i=0; i<10000; i++)
		xy.push_back(x[i]*y[i]);
	double xmedia = CalcolaMedia(x);
	double ymedia = CalcolaMedia(y);
	double xymedia = CalcolaMedia(xy);
	double sx = CalcolaSigma(x);
	double sy = CalcolaSigma(y);
	return (xymedia - xmedia*ymedia)/(sx*sy);
}


//calcola la DEVIAZIONE STANDARD
template <typename T> double CalcolaDev( const vector<T> & v) {
	return sqrt(CalcolaVarianza(v));
};


//Calcola la deviazione standard del campione
template <typename T> double CalcolaDevCampione( const vector<T> & v) {
	double scarto = 0;
	double media = CalcolaMedia(v);
	for (int i=0; i<v.size(); i++)
  	scarto = scarto + pow(v[i] - media, 2);
	return sqrt(scarto/(v.size()-1));
};

template <typename T> double GetMax(const vector<T> &v) {
	double max = 0;
	for(int i=0; i<v.size(); i++)
		if (v[i] > max)
			max = v[i];

	return max;
}

template <typename T> double GetMin(const vector<T> &v) {
	double min = 0;
	for(int i=0; i<v.size(); i++)
		if (v[i] < min)
			min = v[i];

	return min;
}

// =========================================================================
//NUOVE FUNZIONI
//==========================================================================


template <typename T> double error( const vector<T> & AV, const vector<T> & AV2, int n) {
	/*dove 
	n e' il numero di blocchi -1
	AV e' il vettore i cui elementi sono le medie progressive di ciascun blocco
	AV2 e' il vettore i cui elementi sono le medie progressive al quadrato di ciascun blocco	
	*/

	if (n==0)
	  return 0;
	else
	  return sqrt((AV2[n] - pow(AV[n],2))/n);
};

//FUNZIONE CHE SVOLGE LA STATISTICA A BLOCCHI
template <typename T> void BlockStat(int N, const vector<T>& data, vector<T>& ave, vector<T>& ave2, vector<T>& sum_prog, vector<T>& sum2_prog, vector<T>& err_prog) {

	int L = int(data.size()/N);
	ave.clear();
	ave2.clear();
	sum_prog.clear();
	sum2_prog.clear();
	err_prog.clear();

	for(int i=0; i<N; i++) {
    double sum = 0;
    for(int j=0; j<L; j++){ //aggiorno il valore della somma e della somma quadrata
    	double k = j+i*L;
      sum += data.at(k); 
    }
    ave.push_back(sum/L);   
    ave2.push_back(pow((ave.at(i) ),2));
  }

	for(int i=0; i<N; i++) {
    sum_prog.push_back(0);
    sum2_prog.push_back(0);
		for(int j=0; j<i+1; j++) {
    	sum_prog.at(i) += ave.at(j);
    	sum2_prog.at(i) += ave2.at(j);
    }
    
    sum_prog.at(i) /= (i+1); //Cumulative average
    sum2_prog.at(i) /= (i+1); //Cumulative square average
    err_prog.push_back( error(sum_prog, sum2_prog, i) ); //Statistical uncertainty   
  }
};


//STAMPARE I RISULTATI SU FILE
template <typename T>void Print(const vector<T>& sum_prog, const vector<T>& err_prog, int size, const char* NomeFile) {
	
	ofstream out(NomeFile);
	for(int i=0; i<size; i++)
		out << i << " " << sum_prog.at(i) << " " << err_prog.at(i) << endl;
	
	out.close();
};

template <typename T>void PrintV(const vector<T>& vec, const char* NomeFile) {
	
	ofstream out(NomeFile);
	for(long unsigned int i=0; i < vec.size(); i++)
		out << vec.at(i) << endl;
	
	out.close();
};
