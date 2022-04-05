#include <iostream>
#include "Grain.hpp"
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <fstream>

using namespace std;

int N_grain = 10;
double dt = 1e-3;
double xmin = 0., xmax = 1., ymin=0., ymax = 1.;
double g = 9.81;
double dissipation = 0.5;

void initialization (Grain* tab_grain, Grain* tab_grain_copy)
{
	double randNum = 0.;
	double vxmin = 0, vxmax = (xmax-xmin)/(200.*dt);
	double vymin = 0, vymax = (ymax-ymin)/(200.*dt);
	srand(time(NULL));
	
	for (int i = 0; i<N_grain; i++)
		{
			randNum = (double)rand()/(RAND_MAX); //Nombre aléatoire entre 0 et 1
			double x0 = xmin + (xmax-xmin)*randNum;
			
			randNum = (double)rand()/(RAND_MAX);
			double y0 = ymin + (ymax-ymin)*randNum;
			
			Vecteur r0 (x0,y0);
			double rayon = 0.05;
			double rho = 100.;
			
			Grain my_grain(r0,rayon,rho);
			
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


bool est_en_contact(Grain grain1, Grain grain2)
{
	Vecteur pos = grain1.r - grain2.r;
	
	double r1 = grain1.get_rayon();
	double r2 = grain2.get_rayon();
	
	if (pos.norme() < r1+r2)
		return true;
		
	return false;
}


vector<int> contact_list (Grain* tab_grain, Grain & main_grain)

/*
 * @brief Takes a tab of grain, the index of one of them  
 * @return Index of grains in contact with the grain index_grain
 */
 
{
  vector<int> list_contact;

  ;
	for (int i = 0; i < N_grain; i++)
	{
		
	  if (est_en_contact(main_grain, tab_grain[i]) and tab_grain[i].r.get_x() != main_grain.r.get_x())
	    // Does not consider the main_grain in contact_list
	    {
	    list_contact.push_back(i);
	    //cout<< list_contact[0] << " " << list_contact[1]  <<endl;

	    }
	}
	//cout<<"Taille: " <<list_contact.size()<<endl;
	return list_contact;
}



void bord(Grain& my_grain) //Simule un choc élastique avec le bord
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
	 
	}

}

double prod(Vecteur X,Vecteur Y ){
  return (X.get_x()*Y.get_x() + X.get_y()*Y.get_y());

}

void rebond_elastique(Grain &grain1, Grain &grain2)
{
        Vecteur r1 = grain1.r;
	Vecteur v1 = grain1.v;
	double m1 = grain1.get_m();

	Vecteur r2 = grain2.r;
	Vecteur v2 = grain2.v;
	double m2 = grain2.get_m();
	
	Vecteur Pos= grain1.r-grain2.r;
	double Distance= Pos.norme();
	Vecteur v1_f;
	Vecteur v2_f;
	
	if (Distance == 0){
	  v1_f=v1;
	  v2_f=v2;
	  
	}
	else{
	  v1_f = v1 - 2*(m2/(m1+m2))*(prod((v1-v2),(r1-r2))/Distance*Distance)*(r1-r2)  ;
	  v2_f = v2 - 2*(m1/(m1+m2))*(prod((v2-v1),(r2-r1))/Distance*Distance)*(r2-r1);

	}
	grain1.v.set_x(v1_f.get_x());
	grain1.v.set_y(v1_f.get_y());
	grain2.v.set_x(v2_f.get_x());
	grain2.v.set_y(v2_f.get_y());
	//cout<< v1x_f << " et " <<v1x   <<endl;
}

void gravite(Grain& my_grain)
{
	double vy = my_grain.v.get_y();
	my_grain.v.set_y(vy + g*dt);
}

void iteration (Grain* tab_grain, Grain* tab_grain_copy)
{
	for(int i = 0; i<N_grain; i++)
	{
	  bord(tab_grain[i]);
	  // cout<<tab_grain[i].v.get_x()<<endl;
	  vector<int> _contact_list = contact_list(tab_grain,tab_grain[i]);
	  //cout<<_contact_list.size()<<endl;
	  for(int j = 0; j < _contact_list.size(); j++)
	    {
	      if ( i<j){
	      //cout<<"ok3"<<endl;
		  
		rebond_elastique(tab_grain[i], tab_grain[j]);
	      //cout<<"ok2"<<endl;
	      }
	    }
	  
	  //cout<<tab_grain[i].v.get_x()<<endl;
	  gravite(tab_grain[i]);
		
		
	  double x = tab_grain[i].r.get_x();
	  double y = tab_grain[i].r.get_y();
	  double vx = tab_grain[i].v.get_x();
	  double vy = tab_grain[i].v.get_y();
	  
	  
	  //cout<<x << " et " << x + vx*dt <<endl;
	  tab_grain_copy[i].r.set_x(x + vx*dt);
	  tab_grain_copy[i].r.set_y(y + vy*dt);
	  tab_grain_copy[i].v.set_x(vx);
	  tab_grain_copy[i].v.set_y(vy);

	  
	}
	copy_tab(tab_grain,tab_grain_copy,N_grain);
}



int main()
{
	Grain *tab_grain, *tab_grain_copy;
	fstream fich;
	
	tab_grain = (Grain*)malloc(N_grain*sizeof(Grain));
	tab_grain_copy = (Grain*)malloc(N_grain*sizeof(Grain));
	
	initialization(tab_grain,tab_grain_copy);

	Vecteur X(0,2);
	Vecteur Y(3,1);
	cout << prod(X,Y) <<endl; 
	
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
