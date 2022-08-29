#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "GenAlgs.h"
#include "random.h"

using namespace std;

/*========================================		OGGETTI SPECIFICI		=============================================*/

Percorso::Percorso(Random* rnd, int npos) {

	m_rnd = rnd;
	m_npos = npos;
	
	return;
}

Percorso::Percorso(string filename) {
	ifstream in (filename);
	
	m_npos = 0;
	m_pos.clear();
	
	Posizione pos;
	int trash;

	in >> trash;
	
	while (!(in.eof())){
		m_npos ++;		
		in >> pos.x;
		in >> pos.y;
		m_pos.push_back(pos);	
		in >> trash; //controllo che non sia finito il file
	}	
}

Percorso::~Percorso(){}

void Percorso::quadrato(double l) {
	
	m_pos.clear();
	Posizione pos;
	
	for (int i=0; i<m_npos; i++){
		pos.x = l*m_rnd->Rannyu();
		pos.y = l*m_rnd->Rannyu();
		m_pos.push_back(pos);
	}
	return;
}
//ctrl F posizioni
void Percorso::CRF(double r) {
	
	m_pos.clear();
	Posizione pos;
	
	for (int i=0; i<m_npos; i++){
		double ang = m_rnd->Rannyu(0,2*M_PI);
		pos.x = r*cos(ang);
		pos.y = r*sin(ang);
		m_pos.push_back(pos);
	}
}


double Percorso::Distance(int index1, int index2) const {

	return sqrt( pow(m_pos.at(index1).x - m_pos.at(index2).x, 2) + pow( m_pos.at(index1).y - m_pos.at(index2).y , 2) );

}

////QUESTA DA SISTEMARE SECONDO MOMENTO
double Percorso::Misura_Lunghezza_Percorso(Individuo percorso) {
	double cum = 0;
	
	int l = percorso.GetLunghezza();
	
	for (int i=0; i<l ; i++) {
		cum += Distance(percorso.GetGene(i), percorso.GetGene((i+1)%l));
	}
	return cum;
}


void Percorso::Schermo() const {
	for (int i=0; i<m_npos; i++)
		cout << i << " " << m_pos.at(i).x << " " << m_pos.at(i).y << endl;
}

void Percorso::Stampa(string filename) const {
	ofstream out(filename);
	
	for (int i=0; i<m_npos ; i++)
		out << i << " " << m_pos.at(i).x << " " << m_pos.at(i).y << endl;
	
	out.close();
}


/*===================================== ALGORITMI GENETICI ===============================================*/

/*------------------------------------ INDIVIDUO -------------------------------------*/


Individuo::Individuo(Random* rnd, int lunghezza_DNA, vector<int> DNA) {
	m_lDNA = lunghezza_DNA;
	m_rnd = rnd;
	m_DNA = DNA;
	m_fit = 0;
}

Individuo::Individuo(Random* rnd, int lunghezza_DNA) {
	m_lDNA = lunghezza_DNA;
	m_rnd = rnd;
	m_fit = 0;
	for (int i=0; i<m_lDNA; i++)
		m_DNA.push_back(i);
}

Individuo::~Individuo(){}
		
Individuo Individuo::Copia(){
	return Individuo(m_rnd, m_lDNA, m_DNA);
}

int Individuo::GetLunghezza() const{
	return m_lDNA;
}

int Individuo::GetGene(int index) const {
	return m_DNA.at(index);
}

vector<int> Individuo::GetDNA() const{
	return m_DNA;
}

double Individuo::GetFit() const {
	return m_fit;
}

Random* Individuo::GetRND() const {
	return m_rnd;
}

vector<int> Individuo::GetGenes (int index1, int index2) const {
	
	if (index1 > index2)
		swap(index1, index2);
		
	if (index1 > m_lDNA || index2 > m_lDNA) {
		cout << "index out of range 'lunghezza_DNA' " << endl;
		exit(1);
	}
	
	vector<int> Genes (m_DNA.begin()+index1, m_DNA.begin()+index2);
	
	return Genes;
}

void Individuo::operator= (Individuo i) {
	this -> m_lDNA = i.GetLunghezza();
	this -> m_fit = i.GetFit();
	this -> m_rnd = i.GetRND();
	this -> m_DNA = i.GetDNA();	
}

void Individuo::Vincoli(){
	//sposto la prima città all'inizio
	auto it = find (m_DNA.begin(), m_DNA.end(), 0);
	vector<int> cut (it, m_DNA.end());
	m_DNA.erase (it, m_DNA.end());
	m_DNA.insert(m_DNA.begin(), cut.begin(), cut.end());
}

void Individuo::Fit(Percorso* p){
	double fitness = 0;
	for (int i=0; i<m_lDNA-1; i++)
		fitness += p->Distance(m_DNA.at(i), m_DNA.at(i+1));

	m_fit = fitness + p->Distance(m_DNA.at(m_lDNA-1), m_DNA.at(0));
}

void Individuo::Permutazione(int npermutazioni){

	for(int i=0; i<npermutazioni; i++) {
		int x1 = 0; 
		int x2 = 0;	
		while(x1 == x2) {
			x1 = m_rnd->Rannyu(0, m_lDNA);
			x2 = m_rnd->Rannyu(0, m_lDNA);
		}
		swap(m_DNA.at(x1), m_DNA.at(x2));
	}
}

void Individuo::Permutazione(int ext1, int ext2) {
	int delta = ext2-ext1;
	int npermutazioni = m_rnd->Int(1,delta/2);
	
	for(int i=0; i<npermutazioni; i++) {
		int x1 = 0;
		int x2 = 0;		
		
		while(x1 == x2) {
			x1 = m_rnd->Int(ext1, ext2);
			x2 = m_rnd->Int(ext2+1, +delta);
		}
		
		swap(m_DNA.at(x1), m_DNA.at(x2));
	}
};



void Individuo :: Mutazione1() {
	Permutazione(1);
}

void Individuo::Mutazione2() {
	int inizio = 0;
	int fine = 0;
	
	while (inizio == fine){
		inizio = m_rnd->Int(1, m_lDNA/2);
		fine = m_rnd->Int(m_lDNA/2 +1, m_lDNA-1);
	}
	
	if (inizio > fine)
		swap(inizio, fine);

	Permutazione(inizio, fine);

}


void Individuo::Mutazione3() {
	int inizio = 0;
	int fine = 0;
	
	while (inizio == fine){
		inizio = m_rnd->Int(1, m_lDNA-1);
		fine = m_rnd->Int(1, m_lDNA-1);
	}
		
	if (inizio > fine)
		swap(inizio, fine);
		
	if(fine-inizio%2 == 0) { //pari
		while( fine != inizio ) {
			swap(m_DNA.at(inizio), m_DNA.at(fine));
			inizio += 1;
			fine -= 1;
		}
		
	} else { //dispari	
			int i=0;
			while(i<(fine-inizio)/2){
				swap(m_DNA.at(inizio), m_DNA.at(fine));
				inizio += 1;
				fine -= 1;
				i++;
			}
	} 

};

void Individuo::Shift() {
	int primo = 0;
	int ultimo = 0;
	
	while (primo == ultimo){
		primo = m_rnd->Rannyu(0, m_lDNA);
		ultimo = m_rnd->Rannyu(0, m_lDNA);
	}
	
	if (primo > ultimo)
		swap(primo,ultimo);
	
	vector<int> shifting ( m_DNA.begin()+primo, m_DNA.begin()+ultimo+1 );
	reverse( shifting.begin(), shifting.end() ); 
	m_DNA.erase ( m_DNA.begin()+primo, m_DNA.begin()+ultimo+1 );
	m_DNA.insert( m_DNA.begin()+primo, shifting.begin(), shifting.end() );	
}


void Individuo::Schermo(){
	for(int i=0; i<m_lDNA; i++)
		cout << m_DNA.at(i) << " \t";

	cout << " →  " << m_fit << endl;
}

void Individuo::Stampa(string filename){	
	ofstream out (filename);
	for(int i=0; i<m_lDNA; i++)
		out << m_DNA.at(i) << " " << m_DNA[(i+1)%m_lDNA] << endl;
		//li stampo in questo formato in modo da avere sulla prima colonna l'indice di partenza e sulla seconda quello di arrivo
	
	out.close();
}


//rivedere
Individuo Individuo::Riproduzione(Individuo parent2) {
	
	int index = m_rnd -> Int(1,m_lDNA-2);
	vector<int> patrimonio (m_DNA.begin(), m_DNA.begin()+index);
	int j=0;

	while (index < m_lDNA){
		if ( find(patrimonio.begin(), patrimonio.end(), parent2.GetGene(j)) == patrimonio.end() ) {		
			patrimonio.push_back (parent2.GetGene(j));
			index++;
		}
		j++;
	}
	Individuo figlio (m_rnd, m_lDNA, patrimonio);
	return figlio;
}


/*------------------------------------ POPOLAZIONE -------------------------------------*/


Popolazione :: Popolazione(Random* rnd, Percorso* p, int n_individui, int lunghezza_DNA, double select, double mut){
	
	m_rnd = rnd;
	m_p = p;
	m_nindividui = n_individui;
	m_lDNA = lunghezza_DNA;
	m_selettivita = select;
	m_mutabilita = mut;
}

Popolazione :: ~Popolazione(){}


void Popolazione::GeneraPop() {
	
	Individuo primo(m_rnd, m_lDNA);
	primo.Fit(m_p);	
	m_pop.push_back(primo);
	
	for(int i=1;i<m_nindividui;i++) {
		primo.Permutazione(m_lDNA);
		primo.Vincoli();
		primo.Fit(m_p);
		m_pop.push_back(primo);
		
	}
	Order();
}


Individuo Popolazione::GetInd(int index) const {
	return m_pop.at(index);
}

double Popolazione::GetMedia() const {
	return m_media;
}

double Popolazione::GetStdev() const {
	return m_stdev;
}

double Popolazione::GetBest() const {
	return m_best;
}

double Popolazione::GetLunghezza() const {
	return m_lDNA;
}

double Popolazione::GetNindividui() const {
	return m_nindividui;
}

double Popolazione::GetSelect() const {
	return m_selettivita;
}

double Popolazione::GetMut() const {
	return m_mutabilita;
}


void Popolazione :: NewGen(int n_accoppiamenti) {

	Order();
	
	for(int i = 0; i<n_accoppiamenti; i++){
		
		int j = int(ceil(m_nindividui * pow(m_rnd->Rannyu(),m_selettivita)))-1;
		int k = int(ceil(m_nindividui * pow(m_rnd->Rannyu(),m_selettivita)))-1;

		Individuo figlio = m_pop.at(j).Riproduzione(m_pop.at(k));
		
		if (m_rnd->Rannyu() < m_mutabilita){
			
			switch (m_rnd->Int(1,4)){
				case 1: figlio.Mutazione1();
				break;
				case 2: figlio.Mutazione2();
				break;
				case 3: figlio.Mutazione3();
				break;
				case 4: figlio.Shift();
				break;
			}
		}
		figlio.Vincoli();
		figlio.Fit(m_p);
		Sostituisci(figlio);
	}
	
	Order();

	return;
}




void Popolazione::Order() {

	for(int i=0; i<m_nindividui; i++){
		for(int j=i+1; j<m_nindividui; j++) {
			if( m_pop.at(j).GetFit() < m_pop.at(i).GetFit() )
				swap(m_pop.at(i), m_pop.at(j));			
		}
	}	
}

void Popolazione::Sostituisci(Individuo nuovo) {
	
	int x = m_rnd->Int(ceil(m_nindividui/2),m_nindividui-1);	//sostituisce un individuo della metà meno efficiente della popolazione
	m_pop.at(x) = nuovo;
}


void Popolazione::Statistica() {

	Order();
	double half = 0;
	double sum = 0;
	double sum2 = 0;
	
	for(int i=0; i<m_nindividui; i++){
		sum += m_pop.at(i).GetFit();
		sum2 += pow(m_pop.at(i).GetFit(),2);
		
		if(i==int(m_nindividui/2))
			half = sum;
	}
				
	m_media = half / double(int(m_nindividui/2));
	m_best = m_pop.at(0).GetFit();
	m_stdev = sqrt(sum2/double(m_nindividui) - pow(sum/double(m_nindividui),2));
}

void Popolazione::Schermo() {
	cout << "** POPOLAZIONE **" << endl;
	for(int i=0; i<m_nindividui; i++){
		m_pop.at(i).Schermo();
	}
	cout << "*****************" << endl;
}
