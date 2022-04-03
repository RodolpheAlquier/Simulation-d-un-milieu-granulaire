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

void initialization (Grain* tab_grain, Grain* tab_grain_copy)
{
	double randNum = 0.;
	double vxmin = 0, vxmax = (xmax-xmin)/(2.*dt);
	double vymin = 0, vymax = (ymax-ymin)/(2.*dt);
	srand(time(NULL));
	
	for (int i = 0; i<N_grain; i++)
		{
			randNum = (double)rand()/(RAND_MAX + 1); //Nombre aléatoire entre 0 et 1
			double x0 = xmin + (xmax-xmin)/2.*randNum;
			
			randNum = (double)rand()/(RAND_MAX + 1);
			double y0 = ymin + (ymax-ymin)/2.*randNum;
			
			Vecteur r0 (x0,y0);
			double rayon = 1.;
			double rho = 1.;
			Grain my_grain(r0,rayon,rho);
			
			randNum = (double)rand()/(RAND_MAX + 1);
			double vx0 = vxmin + (vxmax-vxmin)/2.*randNum;
			my_grain.v.set_x(vx0);
			
			randNum = (double)rand()/(RAND_MAX + 1);
			double vy0 = vymin + (vymax-vymin)/2.*randNum;
			my_grain.v.set_y(vy0);
			
			
			tab_grain[i] = my_grain;
			tab_grain_copy[i] = tab_grain[i];
			
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
	double x1 = grain1.r.get_x();
	double y1 = grain1.r.get_y();
	double x2 = grain2.r.get_x();
	double y2 = grain2.r.get_y();
	
	double r1 = grain1.get_rayon();
	double r2 = grain2.get_rayon();
	
	if (sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) < r1+r2)
		return true;
		
	return false;
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
		
		if (est_en_contact(main_grain, tab_grain[i]) and i != index_grain) // Does not consider the grain "index_grain" in contact_list
			list_contact.push_back(i);
			
	}
	
	return list_contact;
}


void bord(Grain my_grain) //Simule un choc élastique avec le bord
{
	double x = my_grain.r.get_x();
	double y = my_grain.r.get_y();
	double vx = my_grain.v.get_x();
	double vy = my_grain.v.get_y();
	double rayon = my_grain.get_rayon();
	
	if((xmax-x) <= rayon or (x-xmin) <= rayon)
		my_grain.v.set_x(-vx);  //Choc élastique avec le bord A MODIFIER AVEC LOI HERTZ
	
	if((ymax-y) <= rayon or (y-ymin) <= rayon)
		my_grain.v.set_y(-vy);  //Choc élastique avec le bord
}

void iteration (Grain* tab_grain, Grain* tab_grain_copy)
{
	for(int i = 0; i<N_grain; i++)
	{
		bord(tab_grain[i]);
		
		double x = tab_grain[i].r.get_x();
		double y = tab_grain[i].r.get_y();
		double vx = tab_grain[i].v.get_x();
		double vy = tab_grain[i].v.get_y();
		
		
		
		tab_grain_copy[i].r.set_x(x + vx*dt);
		tab_grain_copy[i].r.set_y(y + vy*dt);
	}
	copy_tab(tab_grain,tab_grain_copy,N_grain);
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








void loi_hertz(Grain grain1, Grain grain2)
{
        double kt=1;
        double kn=1;

	
        double v1x = grain1.v.get_x();
	double v1y = grain1.v.get_y();
	double v1 = sqrt(v1x*v1x + v1y*v1y);
	double m1 = grain1.get_m();
	
	double v2x = grain2.v.get_x();
	double v2y = grain2.v.get_y();
	double v2 = sqrt(v2x*v2x + v2y*v2y);
	double m2 = grain2.get_m();


	double x1 = grain1.r.get_x();
	double y1 = grain1.r.get_y();
	double x2 = grain2.r.get_x();
	double y2 = grain2.r.get_y();
	
	double r1 = grain1.get_rayon();
	double r2 = grain2.get_rayon();

	double d = ((x1-x2)**2+(y1-y2)**2)**0.5;
	double delta_tang = r1+r2-d;
	double delta_norm = 2*(r1**2 - d**2/4 )**0.5; //C'est vraie ssi r1 est equiv r2

	double phi = angle_collision(grain1,grain2);

	double Fx = cos(phi)*kt * delta_tang**1.5 - sin(phi)*kn*delta_norm**1.5 ;
	double Fy = sin(phi)*kt * delta_tang**1.5 + cos(phi)*kn*delta_norm**1.5 ;

	if (v1x*Fx>0){
	  double v1x_f = v1x - Fx*dt/m1;
	}

	if (v1x*Fx<=0){
	  double v1x_f = v1x + Fx*dt/m1;
	}

	if (v2x*Fx<=0){
	  double v2x_f = v2x + Fx*dt/m2;
	}

	if (v2x*Fx>0){
	  double v2x_f = v2x - Fx*dt/m2;
	}

	
	if (v1y*Fy>0){
	  double v1y_f = v1y - Fy*dt/m1;
	}

	if (v1y*Fy<=0){
	  double v1y_f = v1y + Fy*dt/m1;
	}

	if (v2y*Fy<=0){
	  double v2y_f = v2y + Fy*dt/m2;
	}

	if (v2y*Fy>0){
	  double v2y_f = v2y - Fy*dt/m2;
	}



	grain1.v.set_x(v1x_f);
	grain1.v.set_y(v1y_f);
	grain2.v.set_x(v2x_f);
	grain2.v.set_y(v2y_f);


	
}









void rebond_elastique(Grain grain1, Grain grain2)
{
	double v1x = grain1.v.get_x();
	double v1y = grain1.v.get_y();
	double v1 = sqrt(v1x*v1x + v1y*v1y);
	double m1 = grain1.get_m();
	
	double v2x = grain2.v.get_x();
	double v2y = grain2.v.get_y();
	double v2 = sqrt(v2x*v2x + v2y*v2y);
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



int main()
{
	Grain *tab_grain, *tab_grain_copy;
	fstream fich;
	
	tab_grain = (Grain*)malloc(N_grain*sizeof(Grain));
	tab_grain_copy = (Grain*)malloc(N_grain*sizeof(Grain));
	
	initialization(tab_grain,tab_grain_copy);
	
	
	double T = 10*dt;
	fich.open("C:\Users\Antonio Vicente\Desktop\M1\FPT\Bulku\C++\position.csv", ios::out);
	
	for (double t=0; t < T; t = t+dt)
	{
		for(int i = 0; i<N_grain; i++)
		{
			double x = tab_grain[i].r.get_x();
			double y = tab_grain[i].r.get_y();
			double vx = tab_grain[i].v.get_x();
			double vy = tab_grain[i].v.get_y();
			fich << x << " " << y << " " << vx << " " << vy << endl;
		}
		
		iteration(tab_grain,tab_grain_copy);
		
		fich << "###############################################"<< endl;
	}
	fich.close();
	
	return 0;
	
}
