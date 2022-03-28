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
int N_grain=100;                           //Nombre de particules
int constexpr L=100;
int constexpr H=200;                 //Longueur et hauteur de la boîte
double rayon = 10.;
double rho = 1.;



//Definition de la fonction system qui donne l'équation du mvt d'une particule par rapport aux autres



void system (double* q, double* qp, double t)
{
	double x = q[0], y = q[1];
	double v_x = q[2], v_y = q[3];
	
	qp[0] = v_x;
	qp[1] = v_y;
	qp[2] = 1./m * (K_t*pow(delta_t,3./2.)*cos(theta) - K_n*pow(delta_n,3./2.)*sin(theta));
	qp[3] = 1./m * (K_t*pow(delta_t,3./2.)*sin(theta) + K_n*pow(delta_n,3./2.)*cos(theta));
	
	
}



//Definition de la fontion Runge-Kutta
/*
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
*/


	
 //Initialisation du tableau Grains
	
 void initialization (Grain* tab_grain){
   for (int i = 0; i<N_grain; i++){
     Vecteur r0 ( rand() % L , rand() % H );  
     double rayon = 10;
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




void copy_tab(Grain* tab_grain , Grain* tab_grain_syst){
  for(int i; i<N_grain;i++){
    double  X = tab_grain[i].r.get_x();
    double  Y = tab_grain[i].r.get_y();
    double R = tab_grain[i].get_rayon();
    double Rho = tab_grain[i].get_rho();

    Vecteur R0 ( X , Y );  
   
    Grain my_grain(R0,R,Rho);

    tab_grain_syst[i]=my_grain;
  }
}


void le_temps_passe(Grain* tab_grain,double temps, double dt){
  Grain * tab_temp;
  tab_temp = (Grain*)malloc(N_grain*sizeof(Grain));
  for (double t; t<temps; t=t+dt){
    copy_tab(tab_grain,tab_temp);
   
    for(int k=0; k< N_grain; k++){
      tab_temp[k].r.set_x(t);
      if (k==1){
	cout<< tab_temp[k].r.get_x()<<endl;
      }
      //vector<int> _list_pos;
      //_list_pos = contact_list(tab_grain,k);
      // cout<<t<<endl;
      //cout<<_list_pos.size() <<" time: " << t << " particule: " << k  <<endl;
      //changement de position avec la fonction systeme
      // delete [] _list_pos;
    }
    cout<<"Tableau au temps : "<<t<<endl;   


  }


}

void Conditions_limites(Grain* tab_grain,double temps, double dt){
   for(int k; k< N_grain; k++){
     double x_lim = tab_grain[k].get_r().get_x();
     double y_lim = tab_grain[k].get_r().get_y();
     double vx_lim = tab_grain[k].get_v().get_x();
     double vy_lim = tab_grain[k].get_v().get_y();
     if( x_lim > L){
       tab_grain[k].r.set_x(L);
       tab_grain[k].v.set_x(-vx_lim);
     }
      if( x_lim < 0){
       tab_grain[k].r.set_x(0);
       tab_grain[k].v.set_x(-vx_lim);
     }
      if( y_lim > H){
       tab_grain[k].r.set_y(H);
       tab_grain[k].v.set_y(-vy_lim);
     }
      if( x_lim < 0){
       tab_grain[k].r.set_x(0);
       tab_grain[k].v.set_y(-vy_lim);
     }
     
   }	 
       
}







 

 












int main(){

  srand(time(NULL));
  
	
  int iteration = 1000;
  double temps= 1.;
  double dt;
  dt = temps/iteration;


  	
  Grain* tab_grain;
  vector<int> _list;
  
  tab_grain = (Grain*)malloc(N_grain*sizeof(Grain));
  
  initialization(tab_grain);
  
 
  
  _list = contact_list(tab_grain,2);
	 
  cout<<_list.size()<<endl;




  //VERIFICATIONS:
  
  
  //VERIFICATION de la bonne initialisation
  
  /* for(int k=0; k< N_grain; k++){
    cout<< tab_grain[k].r.get_x()<<"  " <<tab_grain[k].r.get_y() <<endl;
   }

  for (int i = 0; i <_list.size(); i++){
    // cout << _list[i]<<"list" << endl;
  } //Changer rayon à 10 pour la verification Ok

  */



 //VERIFICATION que la copie du tableau fonction bien



 /*
Grain * tab_grain_syst ;
tab_grain_syst =(Grain*)malloc(N_grain*sizeof(Grain));
copy_tab(tab_grain,tab_grain_syst);
  
for(int k=0; k< N_grain; k++){
   tab_grain_syst[k].r.set_x(1);
   cout<< tab_grain_syst[k].r.get_x()<<" system  " <<tab_grain[k].r.get_x()<<"original"<<endl;
   }
   // le_temps_passe(tab_grain,30,0.01);

//Ok
*/


  //VERIFICATION de la fontion le temps passe


  //le_temps_passe(tab_grain, temps, dt);

  //Ok
  










  
  //VERIFICATION de l'iteraction
  

  //Iteraction
  double x1,y1,vx1,vy1;
  double x2,y2,vx2,vy2;

  vector<double> q1 (4);
  vector<double> q2 (4);
  vector<double> q_r(2),q_v(2);
  
  


  
  /* for (int t=0; t<iteration; t=t+dt){
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
 }*/


 

  
  return 0;



}



