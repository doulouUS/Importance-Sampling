//
//  random.cpp
//  echantillonnage_preferentiel
//
//  Created by DOUGE Louis on 4/14/16.
//  Copyright © 2016 DOUGE Louis. All rights reserved.
//
#include "random.hpp"
#include "matrice.hpp"

/// Renvoie un nombre aleatoire uniforme
double uniformRandom()
{
    return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. ); // Trouvé sur internet
}

/// Renvoie un nombre aleatoire gaussien
double normalRandom_1()
{
    /// Filtrage des 0 dans Box-Muller
    double u1(0);
    while(u1==0)
    {
        u1=uniformRandom();
    }
    double u2=uniformRandom();
    return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1)); // on utilise qu'une partie de Box-Muller
}

/// Renvoie un vector de x[i] suivant une loi gaussienne

vector<double> normalRandomSimulation(int N)
{
    vector<double> X (N);
    for(int i=0;i<N;i++)
    {
        double x = normalRandom_1(); /// Nombre aleatoire suivant une loi normale
        //cout << " x_" << i << "= " << x << endl;
        X[i] = x;
    }
    
    return X;
}


/// Renvoie trois nombres aléatoires gaussiens indépendants 
vector<double> normalRandom_3()
{
    
    double u1=0;
    double u2=0;
    double u3=0;
    double u4=0;
    
    while (u1==0 || u2==0 || u3==0 || u4==0) {
        u1=uniformRandom();
        u2=uniformRandom();
        u3=uniformRandom();
        u4=uniformRandom();
    }
    
    
    double x1=cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
    double x2=sin(8.*atan(1.)*u1)*sqrt(-2.*log(u2));
      //double x3=sin(8.*atan(1.)*u3)*sqrt(-2.*log(u4));
    double x3=cos(8.*atan(1.)*u4)*sqrt(-2.*log(u3));
    
    vector<double> X (3);
    X[0] = x1;
    X[1] = x2;
    X[2] = x3;
    
    return X;
}

///Renvoie un mouvement brownien de matrice de corrélation définie dans l'énoncé, arrêté au temps T.
vector<double> W_T_gamma (double rho_1, double rho_2, double rho_3, double T, vector<double>& X){
    // A partir d'un echantillon de 3 variables normales independantes X, on calcule le brownien
    // les arguments permettent de former la matrice gamma telle que définie dans l'énoncé
    /*
    /   1       rho_1   rho_2   \
    |   rho_1   1       rho_3   |
    \   rho_2   rho_3   1       /
     
     T est le temps où l'on veut arrêter le mouvement browninen (on n'a pas besoin de toute la trajectoire pour calculer les prix au temps T).
    */
    //matrice de corrélation
    matrice gamma(3,3,1);
    gamma(1,2)=rho_1; gamma(1,3)=rho_2; gamma(2,3)=rho_3;
    gamma(2,1)=rho_1; gamma(3,1)=rho_2; gamma(3,2)=rho_3;
    
    //décomposition de cholesky
    matrice chol=gamma.cholesky_3d();
    
    //Brownien final W_T ~ sqrt(T)*AX où A est la décompo de Cholesky de Gamma
    vector<double> resultat=chol*X;
    
    for (int i=0; i<resultat.size(); i++) {
        resultat[i]*=sqrt(T);
    }
    return resultat;
}

