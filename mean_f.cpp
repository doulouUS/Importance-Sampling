//
//  mean_f.cpp
//  echantillonnage_preferentiel
//
//  Created by DOUGE Louis on 4/14/16.
//  Copyright © 2016 DOUGE Louis. All rights reserved.
//

#include "mean_f.hpp"
#include <vector>
#include <cmath>
//#include "engine.h"
#include "matrice.hpp"
#include "random.hpp"
#

using namespace std;

//------------------------------------------------------------------------------------------------
//
//          CAS UNIDIMENSIONNEL
//
//------------------------------------------------------------------------------------------------

//--------
// Payoff call européen 1D
//--------

double g(double S0, double r, double sigma, double T, double K, double x){
    double result = S0*exp((r-0.5*pow(sigma,2.0))*T+sigma*sqrt(T)*x)-K;
    return  (result<0)?0:result;
};

double f(double S, double r, double sigma, double T, double K, double theta,double x){
    return g(S,r,sigma,T,K,x+theta)*exp(-theta*x-0.5*pow(theta,2.0));
}

//--------
// Espérance de f(theta,X) 1D
//--------

double mean_f (double S0, double r, double sigma, double T, double K, vector<double>& X)
{
    /// Calcul la moyenne de f(theta, X) qui est indépendante de theta
    double mean=0.;
    int N = X.size();
    
    for(int i=0;i<N;i++)
    {
        mean += g(S0, r, sigma, T, K, X[i]);
    }
    mean /= N;
    
    return mean;
};

double mean_f_theta (double S0, double r, double sigma, double T, double K, double theta,vector<double> X)
{
    /// Calcul la moyenne de f(theta, X) (qui est indépendante de theta)
    double mean=0.;
    int N = X.size();
    
    for(int i=0;i<N;i++)
    {
        mean += g(S0, r, sigma, T, K, X[i]+theta)*exp(-theta*X[i]-theta*theta/2);
    }
    mean /= N;
    
    return mean;
}


//--------
// Intervalle de confiance 1D
//--------

double ic_f(double S0, double r, double sigma, double T, double K, double theta, vector<double>& X){
    
    // calcul de la variance
    double var=0;
    int N=X.size();
    
    //calcul de la moyenne
    double mean=mean_f(S0, r, sigma, T, K, X);
    
    for (int i=0; i<N; i++) {
        var+=pow(f(S0, r, sigma, T, K, X[i], theta)-mean, 2.0);
    }
    var= var/(N-1);
    
    // intervalle de confiance
    return sqrt(var/N);
};





//------------------------------------------------------------------------------------------------
//
//          CAS TRIDIMENSIONNEL
//
//------------------------------------------------------------------------------------------------



///generateur de prix
vector<double> prix_3d (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& X){
    
    vector<double> prix_3d(3);
    
    //On genere le brownien arrete au temps T
    vector<double> W=W_T_gamma (rho1, rho2, rho3, T, X);
    
    //On calcule le vecteur prix
    for (int i=0; i<3; i++) {
        prix_3d[i]=S0[i]*exp((r-pow(sigma[i], 2.0)/2)*T+sigma[i]*W[i]);
    }
    
    return prix_3d;
}



/// CALL SUR PANIER

//--------
// Payoff multidimensionnel PANIER
//--------

///payoff pour call sur panier
double g_panier (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, vector<double>& X){
    
    //prix arrete au temps T
    vector<double> prix=prix_3d (S0,r,sigma,T, K, rho1,rho2,rho3, X);
    
    double resultat;
    resultat=lambda[0]*prix[0]+lambda[1]*prix[1]+lambda[2]*prix[2]-K;
    return  (resultat<0)?0:resultat;
}

// f_panier(X,theta)
double f_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, vector<double> lambda,vector<double>& X,vector<double>& theta){
    
    //calcul du vecteur X+theta
    vector<double> XplusTheta (3);
    for (int i=0; i<3; i++) {
        XplusTheta[i]=X[i]+theta[i];
    }
    
    //calcul du produit scalaire
    double prodScal=0.;
    for (int i=0; i<3; i++) {
        prodScal+=X[i]*theta[i];
    }
    
    //calcul de la norme de theta au carre
    double normeCarre=0.;
    for (int i=0; i<3; i++) {
        normeCarre+=pow(theta[i],2.0);
    }
    

    //
    return g_panier(S0, r,sigma,  T, K, rho1, rho2, rho3, lambda, XplusTheta)*exp(-prodScal-0.5*normeCarre);
}

///payoff pour put sur panier (utile pour la dernière question)
double g_panier_put (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, vector<double>& X){
    
    //prix arrete au temps T
    vector<double> prix=prix_3d (S0,r,sigma,T, K, rho1,rho2,rho3, X);
    
    double resultat;
    resultat=K-lambda[0]*prix[0]-lambda[1]*prix[1]-lambda[2]*prix[2];
    return  (resultat<0)?0:resultat;
}

//--------
// Espérance de f_panier(theta,X) PANIER
//--------


double mean_g_panier (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, matrice& X){
    
    /// Calcul la moyenne de g_panier
    double mean=0;
    int N = X.getDim_c();
    
    for(int j=0;j<N;j++)
    {
        vector<double> X_j(3);
        for (int i=0; i<3; i++) {
            X_j[i]=X(i+1,j+1);
        }
        mean += g_panier(S0,r,sigma,T,K,rho1,rho2,rho3, lambda, X_j);
    }
    mean /= N;
    
    return mean;
}

double mean_f_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, matrice& X, vector<double> theta)
{
    /// Calcul la moyenne de f_panier
    double mean=0;
    int N = X.getDim_c();
    
    for(int j=0;j<N;j++)
    {
        vector<double> X_j(3);
        for (int i=0; i<3; i++) {
            X_j[i]=X(i+1,j+1);
        }
        mean += f_panier(S0, r, sigma, T, K, rho1, rho2, rho3,lambda,X_j,theta);
    }
    mean /= N;
    
    return mean;
}

//espérance de g_panier_put
double mean_g_panier_put (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, matrice& X)
{
    /// Calcul la moyenne de g_panier
    double mean=0;
    int N = X.getDim_c();
    
    for(int j=0;j<N;j++)
    {
        vector<double> X_j(3);
        for (int i=0; i<3; i++) {
            X_j[i]=X(i+1,j+1);
        }
        mean += g_panier_put(S0,r,sigma,T,K,rho1,rho2,rho3,lambda, X_j);
    }
    mean /= N;
    
    return mean;
};

//--------
// Intervalle de confiance PANIER
//--------

double ic_f_panier(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, matrice& X, vector<double>& theta){
    
    // calcul de la variance
    double var=0;
    int N=X.getDim_c();
    
    //calcul de la moyenne
    double mean=mean_g_panier(S0, r, sigma, T, K, rho1, rho2, rho3, lambda, X);
    
    //On travaille par colonne de la matrice X
    for (int j=1; j<=N; j++) {
        
        //on extrait le vecteur X_j de la colonne j de X
        vector<double> X_j(3);
        for (int i =1;i<=3 ; i++) {
            X_j[i-1]=X(i,j);
        }
        var+=pow(f_panier(S0, r, sigma, T, K, rho1, rho2, rho3, lambda, X_j,theta)-mean,2.0);
    }
    var= var/(N-1);
    
    // intervalle de confiance
    return sqrt(var/N);
};

//
vector<double> ic_f_panier_put(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, vector<double>& lambda, matrice& X){
    
    vector<double> theta={0,0,0};
    vector<double> resultat (2,0);
    
    // calcul de la variance du put
    double var=0;
    int N=X.getDim_c();
    
    //calcul de la moyenne
    double mean=mean_g_panier_put(S0, r, sigma, T, K, rho1, rho2, rho3, lambda, X);
    
    //On travaille par colonne de la matrice X
    for (int j=1; j<=N; j++) {
        
        //on extrait le vecteur X_j de la colonne j de X
        vector<double> X_j(3);
        for (int i =1;i<=3 ; i++) {
            X_j[i-1]=X(i,j);
        }
        var+=pow(g_panier_put(S0, r, sigma, T, K, rho1, rho2, rho3, lambda, X_j)-mean,2.0);
    }
    var= var/(N-1);
    
    resultat[0]=var;
    
    // calcul de la variance de X
    double var2=0;

    
    //calcul de la moyenne
    double mean2=mean_g_panier(S0, r, sigma, T, K, rho1, rho2, rho3, lambda, X);
    
    //On travaille par colonne de la matrice X
    for (int j=1; j<=N; j++) {
        
        //on extrait le vecteur X_j de la colonne j de X
        vector<double> X_j(3);
        for (int i =1;i<=3 ; i++) {
            X_j[i-1]=X(i,j);
        }
        var2+=pow(f_panier(S0, r, sigma, T, K, rho1, rho2, rho3, lambda, X_j,theta)-mean,2.0);
    }
    var2= var2/(N-1);
    
    resultat[1]=var2;
    
    // intervalle de confiance
    return resultat;
}

///CALL ALTIPLANO


//--------
// Payoff multidimensionnel ALTIPLANO
//--------

///payoff pour call altiplano
double g_altiplano (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3,vector<double>& X){
 
    //prix arrete au temps T
    vector<double> prix=prix_3d (S0,r,sigma,T, K, rho1,rho2,rho3, X);
    
    
    //calcul du terme max
    double max_1=(prix[0]-K<prix[1]-K)?prix[1]-K:(prix[0]-K);//max(S1_T-K,S2_T-K)
    double max_2=prix[2]-K;
    
    double max_altiplano=1.5*((max_1<max_2)?max_2:max_1); // 3/2*max(S1_T-K,S2_T-K,S3_T-K)
    
    //calcul du terme min
    double min_1=(prix[0]-K<prix[1]-K)?prix[0]-K:(prix[1]-K);//min(S1_T-K,S2_T-K)
    double min_2=prix[2]-K;
    
    double min_altiplano=1.5*((min_1<min_2)?min_1:min_2); // 3/2*min(S1_T-K,S2_T-K,S3_T-K)
    
    return 3*K-(prix[0]+prix[1]+prix[2])+min_altiplano+max_altiplano; //3K-(S1_T+S2_T+S3_T)+3/2min(S1_T-K,...)+3/2*max(S1_T-k,...)
    
}

//f_altiplano(X,theta)
double f_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3,vector<double>& X,vector<double>& theta){
    
    //calcul du vecteur X+theta
    vector<double> XplusTheta (3,0);
    for (int i=0; i<3; i++) {
        XplusTheta[i]=X[i]+theta[i];
    }
    
    //calcul du produit scalaire
    double prodScal=0;
    for (int i=0; i<3; i++) {
        prodScal+=X[i]*theta[i];
    }
    
    //calcul de la norme de theta au carre
    double normeCarre=0;
    for (int i=0; i<3; i++) {
        normeCarre+=theta[i]*theta[i];
    }
    
    //g(X+theta)exp(-theta.X-0.5|theta|^2)
    return g_altiplano(S0, r,sigma,  T, K, rho1, rho2, rho3, XplusTheta)*exp(-prodScal-0.5*normeCarre);
}

//--------
// Espérance de f_panier(theta,X) ALTIPLANO
//--------

double mean_g_altiplano (vector<double>& S0, double r, vector<double>& sigma, double T, double K, double rho1, double rho2, double rho3, matrice& X){
    
    /// Calcul la moyenne de g_panier
    double mean=0;
    int N = X.getDim_c();
    
    for(int j=0;j<N;j++)
    {
        vector<double> X_j(3);
        for (int i=0; i<3; i++) {
            X_j[i]=X(i+1,j+1);
        }
        mean += g_altiplano(S0,r,sigma,T,K,rho1,rho2,rho3, X_j);
    }
    mean /= N;
    
    return mean;
}

double mean_f_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, matrice& X, vector<double> theta)
{
    /// Calcul la moyenne de f_panier
    double mean=0;
    int N = X.getDim_c();
    
    for(int j=0;j<N;j++)
    {
        vector<double> X_j(3);
        for (int i=0; i<3; i++) {
            X_j[i]=X(i+1,j+1);
        }
        mean += f_altiplano(S0, r, sigma, T, K, rho1, rho2, rho3,X_j,theta);
    }
    mean /= N;
    
    return mean;
}


double ic_f_altiplano(vector<double> S0, double r, vector<double> sigma, double T, double K, double rho1, double rho2, double rho3, matrice& X, vector<double>& theta){
    
    // calcul de la variance
    double var=0;
    int N=X.getDim_c();
    
    
    //calcul de la moyenne
    double mean=mean_f_altiplano(S0, r, sigma, T, K, rho1, rho2, rho3, X, theta);
    
    
    //On travaille par colonne de la matrice X
    for (int j=1; j<=N; j++) {
        
        //on extrait le vecteur X_j de la colonne j de X
        vector<double> X_j(3);
        for (int i =1;i<=3 ; i++) {
            X_j[i-1]=X(i,j);
        }
        
        var+=pow(f_altiplano(S0, r, sigma, T, K, rho1, rho2, rho3, theta, X_j)-mean,2.0);
        
    }
    var/=(N-1);
    
    // intervalle de confiance
    return sqrt(var/N);
};





