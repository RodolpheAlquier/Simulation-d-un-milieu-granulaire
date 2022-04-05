#include "Grain.hpp"

#include <iostream>
#include <cmath>
using namespace std;

//Constructeurs:Vector
Vecteur::Vecteur(double _x, double _y){
    x=_x;
    y=_y;
    z=0.0;
}

Vecteur::Vecteur(){
    x=0.0;
    y=0.0;
    z=0.0;

}
// Get :POur obtenir la position sans le modifier (RK)
// Set : Pour changer la valeur (move)

double Vecteur::get_x()const{
    return x;
}
double Vecteur::get_y()const{
    return y;
}
double Vecteur::get_z()const{
    return z;
}

void Vecteur::set_x(double abs){
    x=abs;
}
void Vecteur::set_y(double ord){
    y=ord;
}
void Vecteur::set_z(double pro){
    z=pro;
}

//Calcul de la norme d'un vecteur

double Vecteur::norme(){
    return sqrt(x*x+y*y+z*z);
}

//Declaration des operateurs:

Vecteur operator+(const Vecteur& vecteur1,const Vecteur& vecteur2){
    return Vecteur(vecteur1.get_x() + vecteur2.get_x() , vecteur1.get_y() + vecteur2.get_y()  );
}

Vecteur operator-(const Vecteur& vecteur1,const Vecteur& vecteur2){
    return Vecteur(vecteur1.get_x() - vecteur2.get_x() , vecteur1.get_y() - vecteur2.get_y()  );
}

double operator*(const Vecteur& vecteur1,const Vecteur& vecteur2){
    return vecteur1.get_x() * vecteur2.get_x() + vecteur1.get_y() * vecteur2.get_y();
}

Vecteur operator*(double a, const Vecteur& vecteur1){
    return Vecteur(a*vecteur1.get_x() , a* vecteur1.get_y());
}

Vecteur & Vecteur::operator=(const Vecteur & vecteur_)
{
	if(& vecteur_ != this){
	x = vecteur_.get_x();
	y = vecteur_.get_y();
	z = vecteur_.get_z();
	}
	return *this;
	
}


//Constructeurs: Grains

/*Grain::Grain(double x0, double y0, double r0, double rho0){ //Construction d'un point
    x=x0;                           //sans vitesse ni acceleratio
    y=y0;
    v_x=0;
    v_y=0;
    a_x=0;
    a_y=0;
    r=r0;
    rho=rho0;
}
// Get :POur obtenir la position,vitesse,acc sans le modifier (RK)
// Set : Pour changer la valeur (move)
double Grain::get_x()const{
    return x;
}
double Grain::get_y()const{
    return y;
}
void Grain::set_x(double abs){
    x=abs;
}
void Grain::set_y(double ord){
    y=ord;
}
double Grain::get_v_x()const{
    return v_x;
}
double Grain::get_v_y()const{
    return v_y;
}
void Grain::set_v_x(double v_abs){
    v_x=v_abs;
}
void Grain::set_v_y(double v_ord){
    v_y=v_ord;
}
double Grain::get_a_x()const{
    return a_x;
}
double Grain::get_a_y()const{
    return a_y;
}
void Grain::set_a_x(double a_abs){
    a_x=a_abs;
}
void Grain::set_a_y(double a_ord){
    a_y=a_ord;
}
 */




Grain::Grain(Vecteur r0, double rayon0, double rho0){ //Construction d'un point
    r=r0;                           //sans vitesse ni acceleratio
    v=Vecteur();
    a=Vecteur();

    rayon=rayon0;
    rho=rho0;
    m=(4/3.) * M_PI * rho*rayon*rayon*rayon;
}

// Get :Pour obtenir la position,vitesse,acc sans le modifier (RK)
// Set : Pour changer la valeur (move)

double Grain::get_m() const
{
    return m;
}

void Grain::set_m(double masse)
{
    m=masse;
}

double Grain::get_rayon() const
{
	return rayon;
}

double Grain::get_rho() const
{
	return rho;
}

int Grain::get_color() const
{
    return color;
}

void Grain::set_color(int couleur)
{
    color=couleur;
}


bool operator!=(const Grain& grain1,const Grain& grain2){
	
    if (grain1.r.get_x() != grain2.r.get_x() or grain1.r.get_y() != grain2.r.get_y() 
	    or grain1.v.get_x() != grain2.v.get_x() or grain1.v.get_y() != grain2.v.get_y() )
			return true;
			
	return false;
	
}

Grain & Grain::operator=(const Grain & grain_)
{
	if(& grain_ != this){
	r.set_x(grain_.r.get_x()) ;
	r.set_y(grain_.r.get_y()) ;
	v.set_x(grain_.v.get_x()) ;
	v.set_y(grain_.v.get_y()) ;
	a.set_x(grain_.a.get_x()) ;
	a.set_y(grain_.a.get_y()) ;
	
	rayon = grain_.get_rayon();
	m = grain_.get_m();
	rho = grain_.get_rho();
	}
	return *this;
	
}
