#include <iostream>
#include "Grain.hpp"
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <fstream>

using namespace std;

int N_grain = 100;
double dt = 1e-3;
double xmin = 0., xmax = 1., ymin=0., ymax = 1.;
double vFond = 0.;
double g = 10;
double dissipation = 0.999;



bool testCollision (Grain & grain1, Grain & grain2)
{
	double R1 = grain1.get_rayon(), R2 = grain2.get_rayon();
	Vecteur r1 = grain1.r, r2 = grain2.r;
	
	return ((r1-r2).norme()*(r1-r2).norme() <= (R1+R2)*(R1+R2));
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
	  my_grain.r.set_y(ymax-rayon);
	  my_grain.v.set_y(vFond-vy);  //Choc élastique avec le bord A MODIFIER AVEC LOI HERTZ
	  
	}
	
	if(y + rayon > ymax)//Bord inf
	{
	  my_grain.r.set_y(ymax-rayon);
	  my_grain.v.set_y(vFond-vy);  //Choc élastique avec le bord A MODIFIER AVEC LOI HERTZ
	  
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

double positionFond(double A, double omega, double t)
{
	return 1. + A*sin(omega*t);
}

double vitesseFond(double A, double omega, double t)
{
	return A*omega*cos(omega*t);
}

void friction (Grain& grain)
{
	double vx = grain.v.get_x();
	double vy = grain.v.get_y();
	
	grain.v.set_x(dissipation*vx);
	grain.v.set_y(dissipation*vy);
}


void iteration (Grain* tab_grain)
{
	for(int i = 0; i<N_grain; i++)
	  {
	    bord(tab_grain[i]);
	    gravite(tab_grain[i]);
		friction(tab_grain[i]);
	
	    
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
	
	    
	    if (Pos.norme() <= ray1+ray2 +ray1*1e-10 and i<j)
	      {
			  updatePositions(tab_grain[i],tab_grain[j]);
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

	
	double T = 10000*dt;
	fich.open("C:/Users/rodod/OneDrive/Bureau/Informatique/C/position.csv", ios::out);
	
	for (double t=0; t < T; t = t+dt)
	{
		for(int i = 0; i<N_grain; i++)
		{
			double x = tab_grain[i].r.get_x();
			double y = tab_grain[i].r.get_y();
			double rayon = tab_grain[i].get_rayon();
			double vx = tab_grain[i].v.get_x();
			double vy = tab_grain[i].v.get_y();
			fich << x << " " << y << " " << ymax << " " << rayon << endl;
		}
		
		ymax = positionFond(0.1,10.,t);
		vFond = vitesseFond(0.1,10.,t);
		
		
		iteration(tab_grain);
		
		fich << "###############################################"<< endl;
		cout << "Chargement ... " << t/T*100 << " %" << endl;
	}
	fich.close();
	cout<<"Ok0"<<endl;
	return 0;
	
}
