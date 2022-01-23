#include <iostream>
#include "Grain.hpp"
#include <cmath>
#include <vector>
#include <random>
#include <time.h>

using namespace std;

double K_t = 1., K_n = 1. ;          //Coefficients de la loi de Hertz
int n = 4;                           //Taille de vector pour Runge Kutta
double delta_t = 0., delta_n = 0.;   //Parametres du systeme
double m = 0., theta = 0.;           //Masse et angles entre les particules
double d =0.;                        //Distance entre particules
int N=100;                           //Nombre de particules
int constexpr L=100;
int constexpr H=200;                 //Longueur et hauteur de la boîte
double rayon = 1.;


//Definition de la fonction system qui donne l'équation du mvt d'une particule par rapport aux autres


void system (vector<double> q, vector<double> qp, double t)
{
	double x = q[0], y = q[1];
	double v_x = q[2], v_y = q[3];
	
	qp[0] = v_x;
	qp[1] = v_y;

	qp[2] = 0;
	qp[3] = 0;

	for(int j=0; j<=N; j++) {
	  int  x_2= 0; //a definir
	  int y_2=0; // a definir
	  d =pow( ( pow(x-x_2,2)+pow(y-y_2,2) ), 0.5 );
	  if (d<2*rayon){
	    delta_t=2*rayon-d;
	    delta_n=2*sqrt(pow(rayon,2) -pow((d/2),2)) ;
	    theta = atan((y-y_2) /(x-x_2,2));
	    
	    qp[2] += 1./m * (K_t*pow(delta_t,3./2.)*cos(theta) - K_n*pow(delta_n,3./2.)*sin(theta));
	    qp[3] += 1./m * (K_t*pow(delta_t,3./2.)*sin(theta) + K_n*pow(delta_n,3./2.)*cos(theta));
	  }
	
	}
	
}


//Definition de la fontion Runge-Kutta
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


int main(){

  double rayon = 1.;
  double rho = 1.;
  int constexpr L=100;
  int constexpr H=200;                   //Longueur et hauteur de la boîte


  
  N=100;
  srand(time(NULL));
  for (int i=1;i<=N;i++){
    Vecteur r1 ( rand() % L , rand() % H );

    string g = "g" + to_string(i);  
    // Grain g+to_sting(i)(r1,rayon,rho) ?? ;  

  }
  Vecteur r1 (0.,0.);
  Vecteur r2 (1.,0.);
  
  vector<double> x1_list, y1_list;
  vector<double> x2_list, y2_list;
  Grain g1(r1,rayon,rho), g2(r2,rayon,rho);
	
  int iteration = 1000;
  double temps = 1.;
  double dt;
  dt = temps/iteration;
	
  double t;
	
  double x1, y1, vx1, vy1;
  double x2, y2, vx2, vy2;
	
  vector<double> q1 (4), q2 (4);
	
	
	
  t = 0;
  cout<<4<<endl;
	
  return 0;
}
