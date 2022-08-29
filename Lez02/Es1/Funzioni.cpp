#include "Funzioni.h"
#include <iostream>
#include <cmath>

void Coseno::Print() const {
	cout << "f(x) = " << GetA() << "cos(" << GetB() << "x + " << GetC() << ") + " << GetD() << endl;
};

void Parabola::Print () const {
	cout << "f(x) = " << GetA() << "x^2 + " << GetB() << "x + " << GetC() << endl;
};

void Retta::Print() const {
	cout << "f(x) = " << GetA() << "x + " << GetB() << endl;
};
