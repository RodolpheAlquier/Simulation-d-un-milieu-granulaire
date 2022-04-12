#ifndef GRAIN_H
#define GRAIN_H

class Vecteur
{
 private :
  double x;
  double y;
  double z; //Vecteur 2D

 public:
  Vecteur(double _x, double _y);
  Vecteur();

  double get_x() const;
  double get_y() const;
  double get_z() const;
  
  void set_x(double abs);
  void set_y(double ord);
  void set_z(double pro);

  
  double norme();

  friend Vecteur operator+(const Vecteur& vecteur1, const Vecteur& vecteur2);
  friend Vecteur operator-(const Vecteur& vecteur1, const Vecteur& vecteur2);
  friend Vecteur operator*(double a, const Vecteur& vecteur1);
  friend double operator*(const Vecteur& vecteur1, const Vecteur& vecteur2);
  Vecteur & operator=(const Vecteur & vecteur_);

  
  
};



class Grain
{
 public:
  Vecteur r;
  Vecteur v;
  Vecteur a;

  double rayon;
  double rho;
  double m;

 public :

  Grain(Vecteur r0, double rayon0, double rho0);

  double get_rayon() const;
  double get_m() const;
  void set_m(double masse);
  double get_rho() const;
  
  bool operator!=(const Grain& grain){
	
		if (r.get_x() != grain.r.get_x() or r.get_y() != grain.r.get_y() 
			or v.get_x() != grain.v.get_x() or v.get_y() != grain.v.get_y() )
				return true;
			
		return false;
  }
  
  Grain & operator=(const Grain & grain_);
  
  



  
  
};




#endif
