//
//  random.hpp
//  echantillonnage_preferentiel
//
//  Created by DOUGE Louis on 4/14/16.
//  Copyright © 2016 DOUGE Louis. All rights reserved.
//

#ifndef random_hpp
#define random_hpp

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>
#include "matrice.hpp"

using namespace std;

/// Renvoie un nombre aleatoire uniforme
double uniformRandom();

/// Renvoie un nombre aleatoire gaussien
double normalRandom_1();

/// Renvoie un vecteur aleatoire gaussien
vector<double> normalRandomSimulation(int N);

/// Renvoie 3 nombres aleatoires gaussiens (temporaire)
vector<double> normalRandom_3();

///Renvoie la valeur d'une martingale au temps T de matrice de correlation gamma (du type de celle définie dans l'énoncé)
vector<double> W_T_gamma (double rho_1, double rho_2, double rho_3, double T, vector<double>& X);


#endif /* random_hpp */
