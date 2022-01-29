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


void system (vector<double> q, vector<double> qp, double t ,vector<Grain>tab_grain ,int i)
{
	double x = q[0], y = q[1];
	double v_x = q[2], v_y = q[3];

	double x_2,y_2,vx_2,vy_2;
	
	qp[0] = v_x;
	qp[1] = v_y;

	qp[2] = 0;
	qp[3] = 0;

	for(int j=i+1; j<=N; j++) {
	  double x_2= tab_grain[j].get_r().get_x();
	  double y_2=tab_grain[j].get_r().get_y();

	  
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
void rk4 (void (*system)(vector<double>, vector<double>,double,vector<Grain>,int), vector<double> q, double t, double dt,vector<Grain> tab_grain ,int i)
{
	vector<double> p1, p2, p3, p4;
	vector<double> p2_calcul, p3_calcul, p4_calcul;
  
	   system(q,p1,t,tab_grain,i);
	
	for(int i = 0; i<n; i++) 
	{
		p2_calcul[i] = p1[i]*dt/2. + q[i];
	}
	
	system(p2_calcul,p2,t + dt/2.,tab_grain,i);
	
	for(int i = 0; i<n; i++) 
	{
		p3_calcul[i] = p2[i]*dt/2. + q[i];
	}
	
	system(p3_calcul,p3,t + dt/2.,tab_grain,i);
	
	for(int i = 0; i<n; i++) 
	{
		p4_calcul[i] = p3[i]*dt + q[i];
	}
	
	system(p4_calcul,p4,t + dt,tab_grain,i);
	
	
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


  
  int N=100;
  srand(time(NULL));

  //Initialisation du tableau Grains
  
  Grain * tab_grain ;
  tab_grain = new Grain[N];
  for (int i=0;i<N;i++){
    Vecteur r1 ( rand() % L , rand() % H );
    tab_grain[i] = Grain(r1, rayon,rho);  
  }
  
	
  int iteration = 1000;
  double temps= 1.;
  double dt;
  dt = temps/iteration;

  //Verification de la bonne initialisation
  
  for(int k=0; k< N; k++){
    cout<< tab_grain[k].get_r().get_x()<<"  " <<tab_grain[k].get_r().get_y() <<endl;
   }

  //Iteraction
  double x1,y1,vx1,vy1;
  double x2,y2,vx2,vy2;

  vector<double> q1 (4);
  vector<double> q2 (4);
  vector<double> q_r(2),q_v(2);
  
  
  Grain * tab_grain_syst ;
  tab_grain_syst = new Grain[N];

  
  for (int t=0; t<iteration; t=t+dt){
    for (int i=0; i<N ; i++){
      
      x1 = tab_grain[i].get_r().get_x();
      y1 = tab_grain[i].get_r().get_y();
      vx1 = tab_grain[i].get_v().get_x();
      vy1 = tab_grain[i].get_v().get_y();
      
      q1[0] = x1;
      q1[1] = y1;
      q1[2] = vx1;
      q1[3] = vy1;


      q2= rk4 (void (*system)(vector<double>, vector<double>,double,vector<Grain>,int), vector<double> q1, double t, double dt,vector<Grain> tab_grain ,int i);

      q_r[0]= q2[0], q_r[1]=q2[1];
      
      tab_grain_syst[i] = Grain(q_r, rayon,rho);
    }
    
    

  }

  cout<<4<<endl;
	
  return 0;
}
