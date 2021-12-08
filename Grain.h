
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

  Vecteur operator+(const Vecteur& vecteur1, const Vecteur& vecteur2);
  Vecteur operator-(const Vecteur& vecteur1, const Vecteur& vecteur2);
  Vecteur operator*(double a, const Vecteur& vecteur1);
  Vecteur operator*(const Vecteur& vecteur1, const Vecteur& vecteur2);

  
  
};



class Grain
{
 private:
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
  
  
};
