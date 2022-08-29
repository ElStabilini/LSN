#ifndef __Position__
#define __Position__

#include <iostream>
#include <cmath>
#include <string>

using namespace std;

class Position {

public:

  // costruttori
  Position() {
  	m_x = 0; m_y = 0; m_z = 0;
	};
	
  Position(double x, double y, double z) {
		m_x = x;
		m_y = y;
		m_z =z;
	}; 
  // distruttore
  ~Position() {;};

  // metodi
	// Coordinate cartesiane
  double getX() const { return m_x; };       
  double getY() const {	return m_y; };
  double getZ() const { return m_z; };

	// Coordinate sferiche
  double getR() const {return sqrt(pow(m_x, 2) + pow(m_y, 2) + pow(m_z, 2));};       
  double getPhi() const {return atan(m_y/m_x);};
  double getTheta() const {return acos(m_z/getR());};

	//coordinate cilindriche
  double getRho() const {return sqrt(pow(m_x, 2) + pow(m_y, 2)); };     // raggio delle coordinate cilindriche

	//metodo per distanza da un altro punto
  double Distanza(const Position& p) const {
		double dx = p.getX() - m_x;
		double dy = p.getY() - m_y;
		double dz = p.getZ() - m_z;
		return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
	};
	
	Position operator+= (const Position& p) {
		return Position( m_x + p.m_x, m_y + p.m_y, m_z + p.m_z);
	}
	
	Position operator/ (int L) {
		return Position( m_x/L, m_y/L, m_z/L );
	}
	
	/*Position pow(const Position& p, int e) {
		return Position( pow(p.m_x,2), pow(p.m_y,2), pow(p.m_z,2) );
	}*/
	
	Position operator/= (double d) {
		return Position( m_x/d, m_y/d, m_z/d );
	}
	
private:

  double m_x, m_y, m_z;  

};

#endif // __Posizione_h__

