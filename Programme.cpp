#include <iostream>
#include "Grain.hpp"

using namespace std;

int main()
{
	double rayon = 2.;
	double rho = 1.5;
	Vecteur r (1.,1.);
	Grain a (r,rayon,rho);
	
	cout << a.r.get_x()<< endl;
	
	
	return 0;
}