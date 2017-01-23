//
//  newton.cpp
//  echantillonnage_preferentiel
//
//  Created by DOUGE Louis on 4/11/16.
//  Copyright © 2016 DOUGE Louis. All rights reserved.
//

#include "newton.hpp"
#include <math.h>

//--------
// Intervalle de confiance
//--------

//double variance (double theta, int N, vector<double> X){


//    double var=0;
//    for (int i=0; i<N; i++) {
//        var+=pow(X[i]-mean(X), 2.0);
//    }
//    return var/(N-1);
//};


//double f(double S, double r, double sigma, double T, double K, double x,double theta){
//  return g(S,r,sigma,T,K,x+theta)*exp(-theta*x-0.5*pow(theta,2.0));
//}

//----------------------------------------------------------------------
//
//     OPTIMISATION UNIDIMENSIONNELLE
//
//----------------------------------------------------------------------

//--------
// estimateur du gradient de la variance
//--------


// unidimensionnel
double u(double S, double r, double sigma, double T, double K,int N,vector<double>& x,double theta){
    // N nombre de simulation de MC
    // x vecteur des simulations du prix (unidimensionnel ici)
    // theta : parametre reel ici car unidimensionnel
    // S:=S0
    
    double result=0;
    for (int i=0; i<N; i++) {
        result+=(theta-x[i])*pow(g(S,r,sigma,T,K,x[i]),2.0)*exp(-theta*x[i]+0.5*pow(theta,2.0));
    }
    return result/N;
}



//--------
// Gradient de l'estimateur du gradient de la variance (unidimensionnel)
//--------

double grad_u(double S0, double r, double sigma, double T, double K,vector<double>& X, double theta){
    // N nb simulation de MC
    int N = X.size();
    double result(0);
    for (int i=0; i<N; i++) {
        double gg = g(S0, r, sigma, T, K, X[i]);
        result+= pow(gg,2.0)*exp(-theta*X[i]+0.5*pow(theta,2.0))*(1+(theta-X[i])*(theta-X[i]));
    }
    return result/N;
};

//--------
// Programme d'optimisation
//--------


double newton (double S0, double r, double sigma, double T, double K, vector<double>& X){
    // N nombre de simulation MC
    int N=X.size();
    int it_max(1000);
    int nb_it(1);
    double theta = 0;
    double theta_next = 0;
    double grad = grad_u(S0, r, sigma, T, K, X, theta);
    
    while( abs(u(S0, r, sigma, T, K, N, X, theta)) >= 1.0e-10 )
    {
        if(nb_it >= it_max) break;
        
        grad = grad_u(S0, r, sigma, T, K, X, theta);
        theta_next = theta - u(S0, r, sigma, T, K, N, X, theta)/grad;
        
        theta = theta_next;
        //cout << "Iteration : " << nb_it << " Theta = " << theta << endl;
        //cout << "Gradient u = " << grad << endl;
        //cout << "u_N(theta) = " << abs(u(S0, r, sigma, T, K, N, X, theta)) << endl;
        nb_it++;
    }
    
    return theta;
};

//----------------------------------------------------------------------
//
//     OPTIMISATION MULTIDIMENSIONNELLE
//
//----------------------------------------------------------------------

/*
vector<double> u_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, vector<double> theta, double rho1, double rho2, double rho3, vector<double> lambda)

matrice grad_u_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, vector<double> theta, double rho1, double rho2, double rho3, vector<double> lambda)

vector<double> newton_D3(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, double rho1, double rho2, double rho3, vector<double> lambda)
*/

/// PANIER

double dot_product(vector<double> X, vector<double> theta)
{   //calcul du produit scalaire
    if(X.size() != theta.size()) exit(-1);
    double prodScal=0;
    for (int i=0; i<X.size(); i++)
    {
        prodScal+=X[i]*theta[i];
    }
    return prodScal;
}

double carre_norme2(vector<double> theta)
{   //calcul de la norme de theta au carre
    double norme=0;
    for (int i=0; i<3; i++) {
        norme+=theta[i]*theta[i];
    }
    return norme;
}


vector<double> u_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, vector<double> theta, double rho1, double rho2, double rho3, vector<double> lambda)
{
    // N nombre de simulation de MC
    // X matrice des simulations du prix 3xN
    // theta : vecteur de taille 3
    
    vector<double> result(3,0);
    
    vector<double> column_j(3);
    for (int j=1; j<=N; j++) // parcours de toutes les colonnes j de la matrice X
    {
        
        column_j[0]=X(j,1);column_j[1]=X(j,2);column_j[2]=X(j,3);
        
        double g = g_panier(S0, r, sigma, T, K, rho1, rho2, rho3, lambda, column_j); // g_panier de l'échantillon X_i
        double prodScal = dot_product(column_j, theta); //produit scalaire entre le ième échantillon et theta
        double carreNorme = carre_norme2(theta); // norme au carre du vecteur theta
        
        for (int i=0; i<3; i++) // parcours les 3 lignes de l'échantillon
        {
            result[i] += (theta[i]-column_j[i])*g*g*exp(-prodScal+0.5*carreNorme);
        }
    }
    
    for (int i=0; i<3; i++)
    {
        result[i] /= N; // division par N a ne pas oublier
    }
    
    return result;
}


///--------
/// Gradient de l'estimateur de la variance (multidimensionnel)
///--------

matrice grad_u_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, vector<double> theta, double rho1, double rho2, double rho3, vector<double> lambda)
{
    // N nombre de simulation de MC
    // X matrice des simulations du prix 3xN
    // theta : vecteur de taille 3
    
    matrice result(3,3,0);
    
    vector<double> column_j(3);
    for (int j=1; j<=N; j++) // parcours de toutes les colonnes j de la matrice X
    {
        column_j[0]=X(j,1);column_j[1]=X(j,2);column_j[2]=X(j,3);
        
        double g = g_panier(S0, r, sigma, T, K, rho1, rho2, rho3, lambda, column_j); // g_panier de l'échantillon X_i
        double prodScal = dot_product(column_j, theta); //produit scalaire entre le ième échantillon et theta
        double carreNorme = carre_norme2(theta); // norme au carre du vecteur theta
        
        matrice mat_grad_j(3,3,0); // matrice du gradient de (theta-X_j)*exp(-theta.X_j+0.5*theta.theta) pour chaque échantillon Xj
        
        // Remplissage du gradient de u pour l'échantillon X_j
        for(int q=1; q<=3; q++) // parcours des 3 colonnes du gradient
        {
            for(int p=1; p<=3; p++) // parcours des 3 lignes du gradient
            {   // Attention aux décalages d'indices entre les matrices et les vectors
                if(p==q) mat_grad_j(p,q) = 1 + (theta[p-1]-column_j[p-1])*(theta[p-1]-column_j[p-1]); // remplissage spécial sur la diagonale
                else mat_grad_j(p,q) = (theta[p-1]-column_j[p-1])*(theta[q-1]-column_j[q-1]);
            }
            
        }
        // Somme des N matrices de gradient multipliées chacune par g^2*exp(-theta.X_j+0.5theta*theta)
       
        result +=mat_grad_j*(g*g*exp(-prodScal+0.5*carreNorme));
    }
    
    result /= N; // division par N a ne pas oublier
    
    return result;
};



///--------
/// Programme d'optimisation dimension 3
/// Renvoie un vecteur theta star de dimension 3
///--------
vector<double> newton_D3_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, double rho1, double rho2, double rho3, vector<double> lambda)
{
    cout << "Entree dans Newton" << endl;
    int it_max(1000);
    int nb_it=1;
    // Theta_star est ici de dimension 3 (theta_star_1;theta_star_2;theta_star_3)
    vector<double> theta(3,0);
    vector<double> theta_next(3,0);
    
    vector<double> u = u_panier(S0, r, sigma, T, K, N, X, theta, rho1, rho2, rho3, lambda);
    //cout << "Calcul de u fini " << endl;
    
    matrice Grad_u (3,3,0);// = grad_u_panier(S0, r, sigma, T, K, N, X, theta, rho1, rho2, rho3, lambda);
    //cout << "Calcul de Grad fini " << endl;
    
    
    while( sqrt(carre_norme2(u)) >= 1.0e-10 )
    {
        if(nb_it >= it_max) break;
        cout << "Iteration " << nb_it;
        
        Grad_u = grad_u_panier(S0, r, sigma, T, K, N, X, theta, rho1, rho2, rho3, lambda);
        u = u_panier(S0, r, sigma, T, K, N, X, theta, rho1, rho2, rho3, lambda);
        
        //theta_next = theta - u(S0, r, sigma, T, K, N, X, theta)/grad;
        vector<double> d(3); // d = (Grad_u)^-1 * u c'est le déplacement de Newton a chaque itération
        d = (Grad_u.inv_3d())*u;
        
        for(int i=0; i<3; i++)
        {
            theta_next[i] = theta[i] - d[i];
        }
        
        theta = theta_next;
        //cout << "Iteration : " << nb_it << " Theta = " << theta << endl;
        //cout << "Gradient u = " << grad << endl;
        cout << " Norme(u_N(theta)) = " << abs(carre_norme2(u)) << endl;
        nb_it++;
    }
    
    return theta;
};
 

/// ALTIPLANO

vector<double> u_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, vector<double> theta, double rho1, double rho2, double rho3)
{
    // N nombre de simulation de MC
    // X matrice des simulations du prix 3xN
    // theta : vecteur de taille 3
    
    vector<double> result(3,0);
    
    vector<double> column_j(3);
    for (int j=1; j<=N; j++) // parcours de toutes les colonnes j de la matrice X
    {
        
        column_j[0]=X(j,1);column_j[1]=X(j,2);column_j[2]=X(j,3);
        
        double g = g_altiplano(S0, r, sigma, T, K, rho1, rho2, rho3, column_j); // g_panier de l'échantillon X_i
        double prodScal = dot_product(column_j, theta); //produit scalaire entre le ième échantillon et theta
        double carreNorme = carre_norme2(theta); // norme au carre du vecteur theta
        
        for (int i=0; i<3; i++) // parcours les 3 lignes de l'échantillon
        {
            result[i] += (theta[i]-column_j[i])*g*g*exp(-prodScal+0.5*carreNorme);
        }
    }
    
    for (int i=0; i<3; i++)
    {
        result[i] /= N; // division par N a ne pas oublier
    }
    
    return result;
}

///--------
/// Gradient de l'estimateur de la variance (multidimensionnel)
///--------

matrice grad_u_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X, vector<double> theta, double rho1, double rho2, double rho3)
{
    // N nombre de simulation de MC
    // X matrice des simulations du prix 3xN
    // theta : vecteur de taille 3
    
    matrice result(3,3,0);
    
    vector<double> column_j(3);
    for (int j=1; j<=N; j++) // parcours de toutes les colonnes j de la matrice X
    {
        column_j[0]=X(j,1);column_j[1]=X(j,2);column_j[2]=X(j,3);
        
        double g = g_altiplano(S0, r, sigma, T, K, rho1, rho2, rho3, column_j); // g_panier de l'échantillon X_i
        double prodScal = dot_product(column_j, theta); //produit scalaire entre le ième échantillon et theta
        double carreNorme = carre_norme2(theta); // norme au carre du vecteur theta
        
        matrice mat_grad_j(3,3,0); // matrice du gradient de (theta-X_j)*exp(-theta.X_j+0.5*theta.theta) pour chaque échantillon Xj
        
        // Remplissage du gradient de u pour l'échantillon X_j
        for(int q=1; q<=3; q++) // parcours des 3 colonnes du gradient
        {
            for(int p=1; p<=3; p++) // parcours des 3 lignes du gradient
            {   // Attention aux décalages d'indices entre les matrices et les vectors
                if(p==q) mat_grad_j(p,q) = 1 + (theta[p-1]-column_j[p-1])*(theta[p-1]-column_j[p-1]); // remplissage spécial sur la diagonale
                else mat_grad_j(p,q) = (theta[p-1]-column_j[p-1])*(theta[q-1]-column_j[q-1]);
            }
            
        }
        // Somme des N matrices de gradient multipliées chacune par g^2*exp(-theta.X_j+0.5theta*theta)
        
        result +=mat_grad_j*(g*g*exp(-prodScal+0.5*carreNorme));
    }
    
    result /= N; // division par N a ne pas oublier
    
    return result;
};

///--------
/// Programme d'optimisation dimension 3
/// Renvoie un vecteur theta star de dimension 3 ALTIPLANO
///--------
vector<double> newton_D3_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, int N, matrice& X,double rho1, double rho2, double rho3)
{
    cout << "Entree dans Newton" << endl;
    int it_max(1000);
    int nb_it=1;
    // Theta_star est ici de dimension 3 (theta_star_1;theta_star_2;theta_star_3)
    vector<double> theta(3,0);
    vector<double> theta_next(3,0);
    
    vector<double> u = u_altiplano(S0, r, sigma, T, K, N, X, theta, rho1, rho2, rho3);
    //cout << "Calcul de u fini " << endl;
    
    matrice Grad_u (3,3,0);// = grad_u_panier(S0, r, sigma, T, K, N, X, theta, rho1, rho2, rho3, lambda);
    //cout << "Calcul de Grad fini " << endl;
    
    
    while( sqrt(carre_norme2(u)) >= 1.0e-10 )
    {
        if(nb_it >= it_max) break;
        cout << "Iteration " << nb_it;
        
        Grad_u = grad_u_altiplano(S0, r, sigma, T, K, N, X, theta, rho1, rho2, rho3);
        u = u_altiplano(S0, r, sigma, T, K, N, X, theta, rho1, rho2, rho3);
        
        //theta_next = theta - u(S0, r, sigma, T, K, N, X, theta)/grad;
        vector<double> d(3); // d = (Grad_u)^-1 * u c'est le déplacement de Newton a chaque itération
        d = (Grad_u.inv_3d())*u;
        
        for(int i=0; i<3; i++)
        {
            theta_next[i] = theta[i] - d[i];
        }
        
        theta = theta_next;
        //cout << "Iteration : " << nb_it << " Theta = " << theta << endl;
        //cout << "Gradient u = " << grad << endl;
        cout << " Norme(u_N(theta)) = " << abs(carre_norme2(u)) << endl;
        nb_it++;
    }
    
    return theta;
};





