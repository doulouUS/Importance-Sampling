//
//  newton.hpp
//  echantillonnage_preferentiel
//
//  Created by DOUGE Louis on 4/11/16.
//  Copyright © 2016 DOUGE Louis. All rights reserved.
//

#ifndef newton_hpp
#define newton_hpp

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
//#include <Eigen/dense> pour les matrices ?


#include "random.hpp"
#include "mean_f.hpp"
#include "matrice.hpp"

#include <vector>

//--------
// estimateur du gradient de la variance
//--------


// unidimensionnel
double u(double S, double r, double sigma, double T, double K,int N,vector<double>& x,double theta);

//--------
// Gradient de l'estimateur du gradient de la variance (unidimensionnel)
//--------
double grad_u(double S0, double r, double sigma, double T, double K,vector<double>& X, double theta);


//--------
// Programme d'optimisation
//--------


double newton (double S0, double r, double sigma, double T, double K, vector<double>& X);

//----------------------------------------------------------------------
//
//     OPTIMISATION MULTIDIMENSIONNELLE
//
//----------------------------------------------------------------------

/// PANIER

double dot_product(vector<double> X, vector<double> theta);

vector<double> u_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, vector<double> theta, double rho1, double rho2, double rho3, vector<double> lambda);

///--------
/// Gradient de l'estimateur de la variance (multidimensionnel)
///--------

matrice grad_u_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, vector<double> theta, double rho1, double rho2, double rho3, vector<double> lambda);

///--------
/// Programme d'optimisation dimension 3
/// Renvoie un vecteur theta star de dimension 3
///--------
vector<double> newton_D3_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X,double rho1, double rho2, double rho3, vector<double> lambda);


/// ALTIPLANO

vector<double> u_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, vector<double> theta, double rho1, double rho2, double rho3);

matrice grad_u_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, vector<double> theta, double rho1, double rho2, double rho3);

vector<double> newton_D3_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X,double rho1, double rho2, double rho3);


#endif /* newton_hpp */
