# Simulation-d-un-milieu-granulaire


[![forthebadge](http://forthebadge.com/images/badges/built-with-love.svg)](http://forthebadge.com)

Projet C++/Simulation FPT pour l'effet des noix du brésil

### Pré-requis

Pour lancer notre projet vous aurez besoin de:

- Un terminal linux ou un environnement de développement (IDE) en C++
- Un terminal linux , un environnement de développement en python ou accès a jupter notebook

### Installation

Pas nécessaire de réaliser aucune installation


## Démarrage

Pour démarrer notre projet il faut toutes les archives .cpp , .hpp et .ipynb

Tout d'abord il faut executer les programmes .cpp et .hpp avec la commande:

```sh
   g++ LeProgrammeDuTout3.cpp Grain.cpp -o Programme
   ```
Ainsi nous créeons un fichier .csv où se trouvent tous les données de notre simulations. Pour l'affichage nous utilisons jupyter notebook mais d'autres IDE sont possibles pour lancer le code en python.

## Prise en main

- Pour gérer le nombre de grain d'une simulation il faut modifier ```N_grain```. **Attention**, il faudra également modifier  ```N``` dans le code python.
- Pour gérer le pas de temps il faut modifier ```dt```.
- Pour modifier le temps d'une simulation, il faut modifier ```T```. **Attention**, il faudra également modifier  ```steps``` dans le code python.
- Pour gérer les caractéristiques des grains, il faut réaliser les modifications dans la fonction ```initialize```.
- Afin de gérer le nombre de grains bleus, il faut aller dans la fonction ```main``` du code python et modifier la  ligne ```if (ndx <= *Nombre de grains bleus*)```

## Fabriqué avec

- [Jupyter notebook](https://jupyter.org/)
- [Spider](https://www.spyder-ide.org/)
- [CodeLite](https://codelite.org/)


## Auteurs

* **Rodolphe ALQUIER**

* **Antonio VICENTE BECERRIL** 



## License
Licence libre
