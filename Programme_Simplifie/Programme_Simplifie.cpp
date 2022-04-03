#include <iostream>
#include "Grain.hpp"
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <fstream>

using namespace std;

int N_grain = 10;
double dt = 1e-4;
double xmin = 0., xmax = 1., ymin=0., ymax = 1.;

void initialization (Grain* tab_grain, Grain* tab_grain_copy)
{
	double randNum = 0.;
	double vxmin = 0, vxmax = (xmax-xmin)/(200.*dt);
	double vymin = 0, vymax = (ymax-ymin)/(200.*dt);
	srand(time(NULL));
	
	for (int i = 0; i<N_grain; i++)
		{
			randNum = (double)rand()/(RAND_MAX + 1); //Nombre aléatoire entre 0 et 1
			double x0 = xmin + (xmax-xmin)*randNum;
			
			randNum = (double)rand()/(RAND_MAX + 1);
			double y0 = ymin + (ymax-ymin)*randNum;
			
			Vecteur r0 (x0,y0);
			double rayon = 0.1;
			double rho = 1.;
			
			Grain my_grain(r0,rayon,rho);
			
			randNum = (double)rand()/(RAND_MAX + 1);
			double vx0 = vxmin + (vxmax-vxmin)*randNum;
			my_grain.v.set_x(vx0);
			
			randNum = (double)rand()/(RAND_MAX + 1);
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
	Vecteur v = grain1.v - grain2.v;
	
	double r1 = grain1.get_rayon();
	double r2 = grain2.get_rayon();
	
	if (v.norme() < r1+r2)
		return true;
		
	return false;
}


vector<int> contact_list (Grain* tab_grain, Grain main_grain)

/**
 * @brief Takes a tab of grain, the index of one of them  
 * @return Index of grains in contact with the grain index_grain
 */
 
{
	vector<int> list_contact;
	
	
	for (int i = 0; i < N_grain; i++)
	{
		
		if (est_en_contact(main_grain, tab_grain[i]) and tab_grain[i] != main_grain) // Does not consider the main_grain in contact_list
			list_contact.push_back(i);
			
	}
	
	return list_contact;
}


void bord(Grain& my_grain) //Simule un choc élastique avec le bord
{
	double x = my_grain.r.get_x();
	double y = my_grain.r.get_y();
	double vx = my_grain.v.get_x();
	double vy = my_grain.v.get_y();
	double rayon = my_grain.get_rayon();
	
	if((xmax-x) <= rayon and vx > 0)
	{
		my_grain.v.set_x(-vx);  //Choc élastique avec le bord A MODIFIER AVEC LOI HERTZ
	}
	
	if ((x-xmin) <= rayon and vx < 0)
	{
		my_grain.v.set_x(-vx);
	}
	
	if((ymax-y) <= rayon and vy > 0)
	{
		my_grain.v.set_y(-vy);  //Choc élastique avec le bord A MODIFIER AVEC LOI HERTZ
	}
	
	if ((y-ymin) <= rayon and vy < 0)
	{
		my_grain.v.set_y(-vy);
	}

}

double angle_collision(Grain grain1, Grain grain2)
{
	double v1x = grain1.v.get_x();
	double v1y = grain1.v.get_y();
	
	double v2x = grain2.v.get_x();
	double v2y = grain2.v.get_y();
	
	double theta1 = atan(v1y/v1x);
	double theta2 = atan(v2y/v2x);
	
	return theta1 - theta2;
}

void rebond_elastique(Grain grain1, Grain grain2)
{
	double v1x = grain1.v.get_x();
	double v1y = grain1.v.get_y();
	double v1 = grain1.v.norme();
	double m1 = grain1.get_m();
	
	double v2x = grain2.v.get_x();
	double v2y = grain2.v.get_y();
	double v2 = grain2.v.norme();
	double m2 = grain2.get_m();
	
	double theta1 = atan(v1y/v1x);
	double theta2 = atan(v2y/v2x);
	
	double phi = angle_collision(grain1,grain2);
	
	double v1x_f =( (v1*cos(theta1 - phi)*(m1-m2) + 2*m2*v2*cos(theta2-phi))/(m1+m2) ) * cos(phi) - v1*sin(theta1 - phi)*sin(phi);
	double v1y_f =( (v1*cos(theta1 - phi)*(m1-m2) + 2*m2*v2*cos(theta2-phi))/(m1+m2) ) * sin(phi) + v1*sin(theta1 - phi)*cos(phi);
	
	double v2x_f =( (v2*cos(theta2 - phi)*(m2-m1) + 2*m1*v1*cos(theta1-phi))/(m1+m2) ) * cos(phi) - v2*sin(theta2 - phi)*sin(phi);
	double v2y_f =( (v2*cos(theta2 - phi)*(m2-m1) + 2*m1*v1*cos(theta1-phi))/(m1+m2) ) * sin(phi) + v2*sin(theta2 - phi)*cos(phi);
	
	grain1.v.set_x(v1x_f);
	grain1.v.set_y(v1y_f);
	grain2.v.set_x(v2x_f);
	grain2.v.set_y(v2y_f);
}

void iteration (Grain* tab_grain, Grain* tab_grain_copy)
{
	for(int i = 0; i<N_grain; i++)
	{
		bord(tab_grain[i]);
		
		vector<int> _contact_list = contact_list(tab_grain,tab_grain[i]);
			
		for(int j = 0; j < _contact_list.size(); j++)
		{
			rebond_elastique(tab_grain[i], tab_grain[j]);
		}
	
		
		double x = tab_grain[i].r.get_x();
		double y = tab_grain[i].r.get_y();
		double vx = tab_grain[i].v.get_x();
		double vy = tab_grain[i].v.get_y();
		
		
		
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
	
	
	double T = 1000*dt;
	fich.open("C:/Users/rodod/OneDrive/Bureau/Informatique/C/position.csv", ios::out);
	
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
	
	return 0;
	
}
