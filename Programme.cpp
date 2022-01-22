#include <iostream>
#include "Grain.hpp"
#include <cmath>
#include <vector>

using namespace std;

double K_t = 1., K_n = 1. ;
int n = 4;
double delta_t = 0., delta_n = 0.;
double m = 0., theta = 0;

void system (vector<double> q, vector<double> qp, double t)
{
	//double x = q[0], y = q[1];
	double v_x = q[2], v_y = q[3];
	
	qp[0] = v_x;
	qp[1] = v_y;
	qp[2] = 1./m * (K_t*pow(delta_t,3./2.)*cos(theta) - K_n*pow(delta_n,3./2.)*sin(theta));
	qp[3] = 1./m * (K_t*pow(delta_t,3./2.)*sin(theta) + K_n*pow(delta_n,3./2.)*cos(theta));
	
	
}

void rk4 (void (*system)(vector<double>, vector<double>,double), vector<double> q, double t, double dt)
{
	vector<double> p1, p2, p3, p4;
	vector<double> p2_calcul, p3_calcul, p4_calcul;
	
	system(q,p1,t);
	
	for(int i = 0; i<n; i++) 
	{
		p2_calcul[i] = p1[i]*dt/2. + q[i];
	}
	
	system(p2_calcul,p2,t + dt/2.);
	
	for(int i = 0; i<n; i++) 
	{
		p3_calcul[i] = p2[i]*dt/2. + q[i];
	}
	
	system(p3_calcul,p3,t + dt/2.);
	
	for(int i = 0; i<n; i++) 
	{
		p4_calcul[i] = p3[i]*dt + q[i];
	}
	
	system(p4_calcul,p4,t + dt);
	
	
	for(int i = 0; i<n; i++) 
	{
		q[i] += dt/6.*(p1[i] + 2*p2[i] + 2*p3[i] + p4[i]);
	}
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






int main()
{
	Vecteur r1 (0.,0.);
	Vecteur r2 (1.,0.);
	double rayon = 1.;
	double rho = 1.;
	vector<double> x1_list, y1_list;
	vector<double> x2_list, y2_list;
	Grain g1(r1,rayon,rho), g2(r2,rayon,rho);
	
	int iteration = 1000;
	double temps = 1.;
	double dt;
	double t;
	
	double x1, y1, vx1, vy1;
	double x2, y2, vx2, vy2;
	
	vector<double> q1 (4), q2 (4);
	
	dt = temps/iteration;
	
	t = 0;
	
	x1 = get_x(), y1 = r1.get_y();
	x2 = r2.get_x(), y2 = r2.get_y();
	vx1 = r1.get_x(), vy1 = r1.get_y();
	vx2 = r2.get_x(), vy2 = r2.get_y();
	
	x1_list.push_back(x1), y1_list.push_back(y1);
	x2_list.push_back(x2), y2_list.push_back(y2);
	
	q1[0] = x1, q1[1] = y1, q1[3] = vx1, q1[4] = vy1;
	q2[0] = x2, q2[1] = y2, q2[3] = vx2, q2[4] = vy2;
	
	
	for(int i = 1; i<iteration; i++)
	{
		t += dt;
		if(sqrt(x1**2+y**2)+sqrt(x2**2+y2**2) <= 2*rayon)
		{
			rk4(system,q1,t,dt);
			rk4(system,q2,t,dt);
		}
		x1 += vx1*dt, y1 = vy1*dt;
		x2 += vx2*dt, y2 = vy2*dt;
	}
	
	
	
	
	
	
	
	return 0;
}