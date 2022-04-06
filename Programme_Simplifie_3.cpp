#include <iostream>
#include "Grain.hpp"
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <fstream>

using namespace std;

int N_grain = 10;
double dt = 1e-2;
double xmin = 0., xmax = 1., ymin=0., ymax = 1.;
double g = 1;
double dissipation = 0.5;

bool collapse(Grain my_grain,Grain *tab_grain ,int i){
  for(int j=0; j<i;j++){
      Vecteur r1 = tab_grain[j].r;
      Vecteur r2 =  my_grain.r;
      double ray1 =  tab_grain[j].get_rayon();
      double ray2 = my_grain.get_rayon();
      Vecteur Pos = r2-r1;
      
      if (Pos.norme() <= ray1+ray2 )
	{
	  return false;
	}
    }
  return true;
  
}


void initialization (Grain* tab_grain)
{
	double randNum = 0.;
	double vxmin = 0, vxmax = (xmax-xmin);
	double vymin = 0, vymax = (ymax-ymin);
	srand(time(NULL));
	
	for (int i = 0; i<N_grain; i++)
	  {
	    randNum = (double)rand()/(RAND_MAX); //Nombre aléatoire entre 0 et 1
	    double x0 = xmin + (xmax-xmin)*randNum;
			
	    randNum = (double)rand()/(RAND_MAX);
	    double y0 = ymin + (ymax-ymin)*randNum;
			
	    Vecteur r0 (x0,y0);
	    double rayon = 0.05;
	    double rho = 1.;
			
	    Grain my_grain(r0,rayon,rho);
			
	    randNum = (double)rand()/(RAND_MAX);
	    double vx0 = vxmin + (vxmax-vxmin)*randNum;
	    my_grain.v.set_x(vx0);
			
	    randNum = (double)rand()/(RAND_MAX);
	    double vy0 = vymin + (vymax-vymin)*randNum;
	    my_grain.v.set_y(vy0);
			
	    	
	    		
	    while(collapse(my_grain,tab_grain,i) == false ){
	      randNum = (double)rand()/(RAND_MAX); 
	      double x0 = xmin + (xmax-xmin)*randNum;
			
	      randNum = (double)rand()/(RAND_MAX);
	      double y0 = ymin + (ymax-ymin)*randNum;

	      my_grain.r.set_x(x0);
	      my_grain.r.set_y(y0);
	      
	    }
	    tab_grain[i] = my_grain;
	  }
}




void copy_tab(Grain* tab1, Grain* tab2, int N)
/**
 * @brief Copy tab2 in tab1. Works for any matrix of size N
 */
{
	for(int i =0; i < N; i++)
		tab1[i] = tab2[i];
}






void bord(Grain & my_grain) //Simule un choc élastique avec le bord
{
	double x = my_grain.r.get_x();
	double y = my_grain.r.get_y();
	double vx = my_grain.v.get_x();
	double vy = my_grain.v.get_y();
	double rayon = my_grain.get_rayon();
	
	if((xmax-x) <= rayon and vx > 0)//Bord droite
	{
	  my_grain.v.set_x(-vx);  //Choc élastique avec le bord A MODIFIER AVEC LOI HERTZ
	  my_grain.r.set_x(xmax-rayon);
	}
	
	if ((x-xmin) <= rayon and vx < 0)//Bord gauche
	{
	  my_grain.v.set_x(-vx);
	  my_grain.r.set_x(xmin+rayon);
	}
	
	if((ymax-y) <= rayon and vy > 0)//Bord inf
	{
	  
	  my_grain.v.set_y(-vy);  //Choc élastique avec le bord A MODIFIER AVEC LOI HERTZ
	  my_grain.r.set_y(ymax-rayon);
	}
	
	if ((y-ymin) <= rayon and vy < 0)
	{
	  my_grain.v.set_y(-vy);
	  my_grain.r.set_y(ymin+rayon);
	 
	}

}


void gravite(Grain& my_grain)
{
	double vy = my_grain.v.get_y();
	my_grain.v.set_y(vy + g*dt);
}

void iteration (Grain* tab_grain)
{
	for(int i = 0; i<N_grain; i++)
	  {
	    bord(tab_grain[i]);
	    gravite(tab_grain[i]);
	    // cout<<tab_grain[i].v.get_x()<<endl;

	    
	  for(int j = 0; j <N_grain; j++){
	    
	    Vecteur r1 = tab_grain[i].r;
	    Vecteur v1 = tab_grain[i].v;
	    double m1 = tab_grain[i].get_m();
	    double ray1 =tab_grain[i].get_rayon();
	    
	    
	    Vecteur r2 = tab_grain[j].r;
	    Vecteur v2 = tab_grain[j].v;
	    double m2 = tab_grain[j].get_m();
	    double ray2 = tab_grain[j].get_rayon();
	    
	    Vecteur Pos= r1-r2;
	    double Distance= Pos.norme();

	    Vecteur v1_f;
	    Vecteur v2_f;
	
	    
	    if (Pos.norme() <= ray1+ray2 +ray1/1000 and i<j)
	      {
		v1_f = v1 - 2*(m2/(m1+m2))*(( ((v1-v2)*(r1-r2))/(Distance*Distance))*(r1-r2));
		v2_f = v2 - 2*(m1/(m1+m2))*(( ((v2-v1)*(r2-r1))/(Distance*Distance))*(r2-r1));
		
		tab_grain[i].v.set_x(v1_f.get_x());
		tab_grain[i].v.set_y(v1_f.get_y());
		tab_grain[j].v.set_x(v2_f.get_x());
		tab_grain[j].v.set_y(v2_f.get_y());


		

		//cout<<"Ok1"<<endl;
	      }
	  
	  }
	 
	    
	  double x= tab_grain[i].r.get_x();
	  double y=tab_grain[i].r.get_y();
	  double vx= tab_grain[i].v.get_x();
	  double vy= tab_grain[i].v.get_y();

	  tab_grain[i].r.set_x(x + vx*dt);
	  tab_grain[i].r.set_y(y + vy*dt);


	  
	  
	}
}


int main()
{
	Grain *tab_grain, *tab_grain_copy;
	fstream fich;
	
	tab_grain = (Grain*)malloc(N_grain*sizeof(Grain));
	
	initialization(tab_grain);

	Vecteur X(5,2);
	Vecteur Y(3,-1);
	cout << X*Y <<" " << Y.norme() <<endl; 
	
	double T = 1000*dt;
	fich.open("position.csv", ios::out);
	
	for (double t=0; t < T; t = t+dt)
	{
		for(int i = 0; i<N_grain; i++)
		{
			double x = tab_grain[i].r.get_x();
			double y = tab_grain[i].r.get_y();
			double vx = tab_grain[i].v.get_x();
			double vy = tab_grain[i].v.get_y();
			fich << x << " " << y << endl;
		}
		
		iteration(tab_grain);
		
		fich << "###############################################"<< endl;
	}
	fich.close();
	cout<<"Ok0"<<endl;
	return 0;
	
}
