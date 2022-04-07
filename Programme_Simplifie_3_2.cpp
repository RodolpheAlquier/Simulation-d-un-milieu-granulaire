#include <iostream>
#include "Grain.hpp"
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <fstream>

using namespace std;

int N_grain = 100;
double dt = 1e-2;
double xmin = 0., xmax = 1., ymin=0., ymax = 1.;
double g = 10;
double dissipation = 0.5;



bool testCollision (Grain & grain1, Grain & grain2)
{
	double R1 = grain1.get_rayon(), R2 = grain2.get_rayon();
	Vecteur r1 = grain1.r, r2 = grain2.r;
	
	return ((r1-r2).norme()*(r1-r2).norme() <= (R1+R2)*(R1+R2));
}




void initialization (Grain* tab_grain, Grain* tab_grain_copy)
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
	    double rayon = 0.01;
	    double rho = 1.;
			
	    Grain my_grain(r0,rayon,rho);
		
		for(int j = 0; j < i; j++)
		{
			if(testCollision(my_grain,tab_grain[j]))
			{
				randNum = (double)rand()/(RAND_MAX);
				x0 = xmin + (xmax-xmin)*randNum;
			
				randNum = (double)rand()/(RAND_MAX);
				y0 = ymin + (ymax-ymin)*randNum;
				
				my_grain.r.set_x(x0);
				my_grain.r.set_y(y0);
				
				j = -1; //Pour recommencer à 0 sachant qu'il y aura j++ juste après
			}
		}
		
		
			
	    randNum = (double)rand()/(RAND_MAX);
	    double vx0 = vxmin + (vxmax-vxmin)*randNum;
	    my_grain.v.set_x(vx0);
			
	    randNum = (double)rand()/(RAND_MAX);
	    double vy0 = vymin + (vymax-vymin)*randNum;
	    my_grain.v.set_y(vy0);


	    tab_grain[i] = my_grain;
	    tab_grain_copy[i] = my_grain;
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
	my_grain.v.set_y(vy + g*dt); // Car l'axe des y est vers le bas
}

void updatePositions(Grain &grain1, Grain &grain2)
{
	double x1 = grain1.r.get_x(), y1 = grain1.r.get_y(), x2 = grain2.r.get_x(), y2 = grain2.r.get_y();
	double R1 = grain1.get_rayon(), R2 = grain2.get_rayon();
	double distance, overlap;
	
	distance = (grain1.r-grain2.r).norme();
	overlap = distance - R1 - R2;
	
	Vecteur deplacementUnitaire = (1./distance)*(grain2.r - grain1.r);
	
	grain2.r = grain2.r - overlap/2 * deplacementUnitaire;
	grain1.r = grain1.r + overlap/2 * deplacementUnitaire;
	
}

void iteration (Grain* tab_grain, Grain * tab_grain_copy)
{
	for(int i = 0; i<N_grain; i++)
	  {
	    
	    bord(tab_grain[i]); //Pas sur copy
	    gravite(tab_grain[i]); //Pas sur copy

	    tab_grain_copy[i]=tab_grain[i];
	    // cout<<tab_grain[i].v.get_x()<<endl;

	    
	  for(int j = 0; j <N_grain; j++){
	    
	    Vecteur r1 = tab_grain_copy[i].r;
	    Vecteur v1 = tab_grain_copy[i].v;
	    double m1 = tab_grain_copy[i].get_m();
	    double ray1 =tab_grain_copy[i].get_rayon();
	    
	    
	    Vecteur r2 = tab_grain_copy[j].r;
	    Vecteur v2 = tab_grain_copy[j].v;
	    double m2 = tab_grain_copy[j].get_m();
	    double ray2 = tab_grain_copy[j].get_rayon();
	    
	    Vecteur Pos= r1-r2;
	    double Distance= Pos.norme();

	    Vecteur dv1_f;
	    Vecteur dv2_f;
	
	    
	    if (Pos.norme() <= ray1+ray2 +ray1*1e-10 and i<j)
	      {
		updatePositions(tab_grain[i],tab_grain[j]);
		dv1_f =- 2*(m2/(m1+m2))*(( ((v1-v2)*(r1-r2))/(Distance*Distance))*(r1-r2));
		dv2_f =- 2*(m1/(m1+m2))*(( ((v2-v1)*(r2-r1))/(Distance*Distance))*(r2-r1));
		
		double v_i_x = tab_grain[i].v.get_x() ;
		double v_i_y = tab_grain[i].v.get_y() ;
		double v_j_x = tab_grain[j].v.get_x() ;
		double v_j_y = tab_grain[j].v.get_y() ;
		
		tab_grain[i].v.set_x(dissipation* (v_i_x + dv1_f.get_x() ) );
		tab_grain[i].v.set_y(dissipation* (v_i_y + dv1_f.get_y() ) );
		tab_grain[j].v.set_x(dissipation* (v_j_x + dv2_f.get_x() ) );
		tab_grain[j].v.set_y(dissipation* (v_j_y + dv2_f.get_y() ) );
		

		

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
	copy_tab(tab_grain_copy,tab_grain, N_grain);
	
	
}


int main()
{
	Grain *tab_grain, *tab_grain_copy;
	fstream fich;
	
	tab_grain = (Grain*)malloc(N_grain*sizeof(Grain));
	tab_grain_copy = (Grain*)malloc(N_grain*sizeof(Grain));
	
	initialization(tab_grain,tab_grain_copy);

	
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
		
		iteration(tab_grain,tab_grain_copy);
		
		fich << "###############################################"<< endl;
	}
	fich.close();
	cout<<"Ok0"<<endl;
	return 0;
	
}
