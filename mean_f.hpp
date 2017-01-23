//
//  mean_f.hpp
//  echantillonnage_preferentiel
//
//  Created by DOUGE Louis on 4/14/16.
//  Copyright © 2016 DOUGE Louis. All rights reserved.
//

#ifndef mean_f_hpp
#define mean_f_hpp

#include <stdio.h>
#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>
#include "matrice.hpp"

using namespace std;


// Contient tous les payoffs (européens, panier, altiplano) fonction de theta ou non
// Contient aussi les fonctions renvoyant des estimations de l'espérance (moyenne) pour ces payoffs



//------------------------------------------------------------------------------------------------
//
//          CAS UNIDIMENSIONNEL
//
//------------------------------------------------------------------------------------------------

double g(double S0, double r, double sigma, double T, double K, double x);

double f(double S, double r, double sigma, double T, double K, double theta,double x);

double mean_f (double S0, double r, double sigma, double T, double K, vector<double>& X);

double mean_f_theta (double S0, double r, double sigma, double T, double K, double theta,vector<double> X);

//double mean(vector<double>& X);

//intervalle de confiance

double ic_f(double S0, double r, double sigma, double T, double K, double theta, vector<double>& X);

//------------------------------------------------------------------------------------------------
//
//          CAS TRIDIMENSIONNEL
//
//------------------------------------------------------------------------------------------------

/// CALL SUR PANIER

//generateur de prix selon le brownien décrit dans l'énoncé
vector<double> prix_3d (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& X);

double g_panier (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, vector<double>& X);

double f_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, vector<double> lambda,vector<double>& X,vector<double>& theta);

//payoff d'un put sur panier
double g_panier_put (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, vector<double>& X);

double mean_g_panier (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, matrice& X);

double mean_f_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, matrice& X, vector<double> theta);

double mean_g_panier_put (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, matrice& X);

//intervalle de confiance
double ic_f_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, matrice& X,vector<double>& theta);

vector<double> ic_f_panier_put(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, matrice& X);


/// CALL ALTIPLANO

double g_altiplano (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3,vector<double>& X);

double f_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3,vector<double>& X,vector<double>& theta);

double mean_g_altiplano (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, matrice& X);

double mean_f_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, matrice& X, vector<double> theta);

//intervalle de confiance
double ic_f_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, matrice& X,vector<double>& theta);

#endif /* mean_f_hpp */
