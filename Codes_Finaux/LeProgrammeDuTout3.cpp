#include <iostream>
#include "Grain.hpp"
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <fstream>

using namespace std;


//Variables générales:

int N_grain = 1000;
double dt = 1e-3;
double xmin = 0., xmax = 1, ymin=0., ymax = 4.;
double g = 10;
double dissipation = 0.99;
double A=2;
double omega = 10;
double pi= M_PI;
double T = 10000*dt;
double phi=0;


//Savoir si deux particules se touchent ou pas:

bool testCollision (Grain & grain1, Grain & grain2)
{
	double R1 = grain1.get_rayon(), R2 = grain2.get_rayon();
	Vecteur r1 = grain1.r, r2 = grain2.r;

	return ((r1-r2).norme()*(r1-r2).norme() <= (R1+R2)*(R1+R2));
}



//Initialisation du tableau Grain
void initialization (Grain* tab_grain)
{
	double randNum = 0.;
	double vxmin = 0, vxmax = 0;
	double vymin = 0, vymax = (ymax-ymin);
	srand(time(NULL));
	
	for (int i = 0; i<N_grain; i++)
	  {
	    randNum = (double)rand()/(RAND_MAX); //Nombre aléatoire entre 0 et 1
	    double x0 = xmin + (xmax-xmin)*randNum; 

	    randNum = (double)rand()/(RAND_MAX);
	    double y0 = (ymax-ymin) + (ymax-ymin)*randNum;

	    Vecteur r0 (x0,y0);
	    double rayon = 0.02;
	    double rho = 1.;

	    //cout<< "Je passe " <<endl;

	    // Création d'une bille de grand volume
	    /* if (i==0){
	      rayon = 0.2;
	      x0= xmax*0.5;
	      y0= ymax-rayon;
	      }*/

	    Grain my_grain(r0,rayon,rho);


	    //On vérifie qu'il n'y a pas billes que se chévauchent
	    
	    for(int j = 0; j < i; j++)
	      {
		if(testCollision(my_grain,tab_grain[j]))
		  {
		    randNum = (double)rand()/(RAND_MAX); //On recalcul un autre nombre aléatoire
		    x0 = xmin + (xmax-xmin)*randNum;

		    randNum = (double)rand()/(RAND_MAX);
		    y0 =(ymax-ymin) + (ymax-ymin)*randNum;

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



void bord(Grain & my_grain) //Simule un choc élastique avec le bord
{
	double x = my_grain.r.get_x();
	double y = my_grain.r.get_y();
	double vx = my_grain.v.get_x();
	double vy = my_grain.v.get_y();
	double rayon = my_grain.get_rayon();

	if((xmax-x) <= rayon and vx > 0){         //Bord droite
	  my_grain.v.set_x(-vx);                  //Choc élastique avec le bord
	  my_grain.r.set_x(xmax-rayon);
	  my_grain.v.set_y(vy*0.5);
	}

	if ((x-xmin) <= rayon and vx < 0){     //Bord gauche
	  my_grain.v.set_x(-vx);
	  my_grain.v.set_y(vy*0.5);            //Frottement normales
	  my_grain.r.set_x(xmin+rayon);
	}

	if((ymax-y) <= rayon and vy > 0){      //Bord inf
	  my_grain.r.set_y(ymax-rayon);
	  my_grain.v.set_x(vx*0.5);
	  my_grain.v.set_y(-vy); 
	}

	if((y+rayon) > ymax ){      //Bord inf
	  my_grain.r.set_y(ymax-rayon);
	  my_grain.v.set_x(vx*0.5);
	  my_grain.v.set_y(-vy); 
	}

	if ((y-ymin) <= rayon and vy < 0){       //Bord sup
	  my_grain.v.set_y(-vy);
	  my_grain.r.set_y(ymin+rayon);
	}
}




void gravite(Grain& my_grain){
  
	double vy = my_grain.v.get_y();
	my_grain.v.set_y(vy + g*dt); // Car l'axe des y est vers le bas
}




void updatePositions(Grain &grain1, Grain &grain2){
  
	double x1 = grain1.r.get_x(), y1 = grain1.r.get_y(), x2 = grain2.r.get_x(), y2 = grain2.r.get_y();
	double R1 = grain1.get_rayon(), R2 = grain2.get_rayon();
	double distance, overlap;

	distance = (grain1.r-grain2.r).norme();
	overlap = distance - R1 - R2; //On calcule l'overlap entre les deux billes lorsqu'elles se chévauchent

	Vecteur deplacementUnitaire = (1./distance)*(grain2.r - grain1.r);

	grain2.r = grain2.r - overlap/2 * deplacementUnitaire;
	grain1.r = grain1.r + overlap/2 * deplacementUnitaire;
}



void acc_Boite(Grain & my_grain, double t){
  
  double vy = my_grain.v.get_y();
  double ay = A*omega*omega*cos(omega*(t-phi));
  my_grain.v.set_y(vy - ay*dt);
  
}



void friction (Grain& grain){
  
  double vx = grain.v.get_x();
  double vy = grain.v.get_y();
  grain.v.set_x(dissipation*vx);
  grain.v.set_y(dissipation*vy);
  
}



void iteration (Grain* tab_grain,double t)
{
	for(int i = 0; i<N_grain; i++)
	  {
	    bord(tab_grain[i]);

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

	    Vecteur dv1_f;
	    Vecteur dv2_f;


	    if (testCollision(tab_grain[i],tab_grain[j])  and i<j)
	      {
		updatePositions(tab_grain[i],tab_grain[j]);
		dv1_f =dv1_f - 2*(m2/(m1+m2))*(( ((v1-v2)*(r1-r2))/(Distance*Distance))*(r1-r2));
		dv2_f =dv2_f - 2*(m1/(m1+m2))*(( ((v2-v1)*(r2-r1))/(Distance*Distance))*(r2-r1));
	      }

	    tab_grain[i].v.set_x(v1.get_x()+dv1_f.get_x());
	    tab_grain[i].v.set_y(v1.get_y()+dv1_f.get_y());
	    tab_grain[j].v.set_x(v2.get_x()+dv2_f.get_x());
	    tab_grain[j].v.set_y(v2.get_y()+dv2_f.get_y());

	    
	  }

	  gravite(tab_grain[i]);
	  friction(tab_grain[i]);
	  //acc_Boite(tab_grain[i],t,phi);

	  
	    //Mouvement de la boîte
	   
	  if (0.5*T/10+pi/omega  >=t and t >=0.5*T/10){
	    phi=0.5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }
	    
	  if (1.5*T/10+pi/omega  >=t and t >=1.5*T/10){
	    phi=1.5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  if (2.5*T/10+pi/omega  >=t and t >=2.5*T/10){
	    phi=2.5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  if (3.5*T/10+pi/omega  >=t and t >=3.5*T/10){
	    phi=3.5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  if (4.5*T/10+pi/omega  >=t and t >=4.5*T/10){
	    phi=4.5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  if (5.5*T/10+pi/omega  >=t and t >=5.5*T/10){
	    phi=5.5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  if (6.5*T/10+pi/omega  >=t and t >=6.5*T/10){
	    phi=6.5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  if (7.5*T/10+pi/omega  >=t and t >=7.5*T/10){
	    phi=7.5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  if (8.5*T/10+pi/omega  >=t and t >=8.5*T/10){
	    phi=8.5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  if (9.5*T/10+pi/omega  >=t and t >=9.5*T/10){
	    phi=9.5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }
		 
	  if (T/10+pi/omega  >=t and t >=T/10){
	    phi=T/10;
	    acc_Boite(tab_grain[i],t);
	  }
	     
	  if (2*T/10+pi/omega  >=t and t >=2*T/10){
	    phi=2*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  if (3*T/10+ pi/omega  >=t and t >=3*T/10){
	    phi=3*T/10;
	    acc_Boite(tab_grain[i],t);
	  }
	    
	  if (4*T/10+pi/omega  >=t and t >=4*T/10){
	    phi=4*T/10;
	    acc_Boite(tab_grain[i],t);
	  }
	
	  if (5*T/10+pi/omega  >=t and t >=5*T/10){
	    phi=5*T/10;
	    acc_Boite(tab_grain[i],t);
	  }
	    
	  if (6*T/10+pi/omega  >=t and t >=6*T/10){
	    phi=6*T/10;
	    acc_Boite(tab_grain[i],t);
	  }
	     
	  if (7*T/10+pi/omega  >=t and t>=7*T/10){
	    phi=7*T/10;
	    acc_Boite(tab_grain[i],t);
	  }
	    	    
	  if (8*T/10+pi/omega  >=t and t>=8*T/10){
	    phi=8*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  if (9*T/10+pi/omega  >=t and t >=9*T/10){
	    phi=9*T/10;
	    acc_Boite(tab_grain[i],t);
	  }

	  
	  double x= tab_grain[i].r.get_x();
	  double y= tab_grain[i].r.get_y();
	  double vx= tab_grain[i].v.get_x();
	  double vy= tab_grain[i].v.get_y();

	  tab_grain[i].r.set_x(x + vx*dt);
	  tab_grain[i].r.set_y(y + vy*dt);
	}
}


int main()
{


	Grain *tab_grain;
	fstream fich;

	tab_grain = (Grain*)malloc(N_grain*sizeof(Grain));
	initialization(tab_grain);


	fich.open("position.csv", ios::out);

	for (double t=0; t < T; t = t+dt){
		for(int i = 0; i<N_grain; i++){
		  
			double x = tab_grain[i].r.get_x();
			double y = tab_grain[i].r.get_y();
			double rayon = tab_grain[i].get_rayon();
			double vx = tab_grain[i].v.get_x();
			double vy = tab_grain[i].v.get_y();
			double rho = tab_grain[i].get_rho();
			fich << x << " " << y << " " << ymax << " " << rayon <<" "<< rho <<endl;
		}

		iteration(tab_grain,t);

		fich << "###############################################"<< endl;
		cout << "Chargement ... " << t/T*100 << " %" << endl;
	}
	fich.close();
	cout<<"Ok0"<<endl;
	return 0;

}
