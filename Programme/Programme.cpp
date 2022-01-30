#include <iostream>
#include "Grain.hpp"
#include <cmath>
#include <vector>

using namespace std;

double K_t = 1., K_n = 1. ;
int n = 4, N_grain = 100;
double delta_t = 0.1, delta_n = 0.1;
double m = 1.;

void system (double* q, double* qp, double t, double theta)
{
	double x = q[0], y = q[1];
	double v_x = q[2], v_y = q[3];
	
	qp[0] = v_x;
	qp[1] = v_y;
	qp[2] = 1./m * (K_t*pow(delta_t,3./2.)*cos(theta) - K_n*pow(delta_n,3./2.)*sin(theta));
	qp[3] = 1./m * (K_t*pow(delta_t,3./2.)*sin(theta) + K_n*pow(delta_n,3./2.)*cos(theta));
	
	
}

void rk4 (void (*system)(double*, double*,double, double), double* q, double t, double dt)
{
	double *p1, *p2, *p3, *p4;
	vector<double> p2_calcul, p3_calcul, p4_calcul;
	
	p1 = (double*)malloc(4*sizeof(double));
	p2 = (double*)malloc(4*sizeof(double));
	p3 = (double*)malloc(4*sizeof(double));
	p4 = (double*)malloc(4*sizeof(double));
	
	system(q,p1,t);
	
	for(int i = 0; i<4; i++) 
	{
		p2_calcul[i] = p1[i]*dt/2. + q[i];
	}
	
	system(p2_calcul,p2,t + dt/2.);
	
	for(int i = 0; i<4; i++) 
	{
		p3_calcul[i] = p2[i]*dt/2. + q[i];
	}
	
	system(p3_calcul,p3,t + dt/2.);
	
	for(int i = 0; i<4; i++) 
	{
		p4_calcul[i] = p3[i]*dt + q[i];
	}
	
	system(p4_calcul,p4,t + dt);
	
	
	for(int i = 0; i<4; i++) 
	{
		q[i] += dt/6.*(p1[i] + 2*p2[i] + 2*p3[i] + p4[i]);
	}
	
	free(p1);
	free(p2);
	free(p3);
	free(p4);
}

/*void system (double* q, double *qp, double t, double theta, double delta_t, double delta_n, double m)
{
	double x = q[0], y = q[1];
	double v_x = q[2], v_y = q[3];
	
	qp[0] = v_x;
	qp[1] = v_y;
	qp[2] = 1./m * (K_t*pow(delta_t,3./2.)*cos(theta) - K_n*pow(delta_n,3./2.)*sin(theta));
	qp[3] = 1./m * (K_t*pow(delta_t,3./2.)*sin(theta) + K_n*pow(delta_n,3./2.)*cos(theta));
	
	
}

void rk4 (void (*system)(double, double,double,double,double,double,double), double* q, double t, double dt)
{
	double *p1, *p2, *p3, *p4;
	double *p2_calcul, *p3_calcul, *p4_calcul;
	
	p1 = (double*)malloc(n*sizeof(double))
	p2 = (double*)malloc(n*sizeof(double))
	p3 = (double*)malloc(n*sizeof(double))
	p4 = (double*)malloc(n*sizeof(double))
	
	for(int i = 0, i<n, i++) 
		{
		
		}
	
	system(q,p1,t,theta,delta_t,delta_n,m,n)
}*/

void initialization (Grain* tab_grain)
{
	for (int i = 0; i<N_grain; i++)
		{
			Vecteur r0 (i,0);   //Lignes de grains A REMPLACER PAR X ET Y ALEATOIRE
			double rayon = 1.;
			double rho = 1.;
			Grain my_grain(r0,rayon,rho);
			
			tab_grain[i] = my_grain;
			
		}
	
	
}



vector<int> contact_list (Grain* tab_grain, int index_grain)

/**
 * @brief Takes a tab of grain, the index of one of them  
 * @return Index of grains in contact with the grain index_grain
 */
 
{
	vector<int> list_contact;
	
	Grain main_grain = tab_grain[index_grain];
	double x = main_grain.r.get_x();
	double y = main_grain.r.get_y();
	double r = main_grain.get_rayon();
	
	for (int i = 0; i < N_grain; i++)
	{
		double X = tab_grain[i].r.get_x();
		double Y = tab_grain[i].r.get_y();
		double R = tab_grain[i].get_rayon();
		
		if (sqrt((X-x)*(X-x)+(Y-y)*(Y-y)) < R+r and i != index_grain) // Does not consider the grain "index_grain" in contact_list
			list_contact.push_back(i);
			
	}
	
	return list_contact;
}


void contact (Grain* tab_grain, Grain* tab_grain_copy, int index_grain)
{
	vector<int> list_contact = contact_list(tab_grain, index_grain); 
}






int main()
{

	
	Grain* tab_grain;
	vector<int> _list;
	
	tab_grain = (Grain*)malloc(N_grain*sizeof(Grain));
	
	initialization(tab_grain);
	
	for(int i = 0; i < N_grain; i++)
	{
		//cout << tab_grain[i].r.get_x() << endl;
	}
	
	 _list = contact_list(tab_grain,2);
	 
	 for (int i = 0; i <_list.size(); i++)
		 cout << _list[i] << endl;
	
	
	return 0;
}
