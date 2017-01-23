//
//  main.cpp
//  echantillonnage_preferentiel
//
//  Created by DOUGE Louis on 4/7/16.
//  Copyright © 2016 DOUGE Louis. All rights reserved.
//

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <fstream>


#include "random.hpp"
#include "mean_f.hpp"
#include "matrice.hpp"
#include "newton.hpp"

#include <vector>
#include <ctime> //connaître le temps de nos fonctions

using namespace std;

int main() {
    
    //----------------Valeurs utilisés pour les tests----------------------

    int N=1e5; //nb simulation MC
    vector<double> S0={1,1,1};
    double r=0.01;
    double T=1;
    double K=1.25;
    vector<double> sigma={0.25,0.28,0.3};
    
    //matrice de corrélation
    double rho1=0.5;
    double rho2=0.5;
    double rho3=0.5;
    
    //lambda
    vector<double> lambda={1./3.,1./3.,1./3.};
    
    //theta
    double theta=1.35946;//theta optimal selon newton en 1D
    vector<double> theta_3d={0,0,0};
   
    //----------------------------------------------------------------------
    
 
    
    //-----------------------------------------
    // UNIDIMENSIONNEL
    //-----------------------------------------

    /*
    vector<double> simulationsNormal =normalRandomSimulation(N);
    
    vector<double> echantNormal= normalRandom_3();
    vector<double> echantNormal2= normalRandom_3();
    */
    //double xx= g(S, r, sigma, T, K, x[1]);
    //cout << "xx=" << xx << endl;
    
    //double gradd=grad_u(S, r, sigma, T, K,theta,xx);
    //cout << "gradd=" << gradd << endl;
    
    //double estim_grad_var=u(S,r,sigma,T,K,n,x,theta);
    //cout << "estim_grad_var=" << estim_grad_var <<    endl;
    
    //double theta=newton(S0[0], r, sigma[0], T, K,simulationsNormal);
    //cout << "theta=" << theta << endl;
    
    
     //double ff=f(S0[0],r,sigma[0],T,K,theta,simulationsNormal[0]);
     //cout<<ff<< endl;
    
    
     //Genere le prix et son intervalle de confiance
    /*
     double prix=exp(-r*T)*mean_f(S0[0], r, sigma[0], T, K, simulationsNormal);
     cout << "Prix= "<< prix << endl;
     
     
    double var=ic_f(S0[0],r,sigma[0],T,K,theta, simulationsNormal);
    cout << "racine de la variance sur N= " << var << endl ;
    
    */
    
    /*
    matrice gamma(3,3,1);
    gamma(1,2)=rho1; gamma(1,3)=rho2; gamma(2,3)=rho3;
    gamma(2,1)=2; gamma(3,1)=88; gamma(3,2)=23;
    
    cout << gamma << endl;
    
    cout << gamma.transpose_3d();
    
    matrice racine = gamma.cholesky_3d();
    vector<double> X(3);
    
    for (int i=0; i<3; i++) {
        X[i]=racine(3,i+1);
    }
    
    vector<double> gamma_1c=racine*X;
    cout << " " << gamma_1c[0] << " " <<  gamma_1c[1] <<" " << gamma_1c[2] << endl;
    
    */
    
    //-----------------------------------------
    // MULTIDIMENSIONNEL
    //-----------------------------------------

    
    //-----------------------------------------
    //creation matrice des simulations 3xN
    //-----------------------------------------
    clock_t temps;
    matrice simulationNormale(3,N,0);
    for (int j=1; j<=N; j++) {
        //pour chaque colonne on introduit un vecteur simulé par normalRandom_3()
        vector<double> simulationNormaleColonne=normalRandom_3();
        for (int i=0; i<3; i++) {
            simulationNormale(i+1,j)=simulationNormaleColonne[i];
        }
    }
    //cout<<simulationNormale;
    /*
    ofstream myfile1;
    myfile1.open (".\valeurs_normales_1.txt");
    myfile1 << "Premier tirage" << endl;
    for (int i=0; i<N; i++) {
        myfile1 << simulationNormale(1,i)<<endl;
    }

    myfile1.close();
     */
     
    /*
    ofstream myfile2;
    myfile2.open ("valeurs_normales_2.txt");
    myfile2 << "Second tirage tirage" << endl;
    for (int i=0; i<N; i++) {
        myfile2 << simulationNormale(2,i)<<endl;
    }
   
    myfile2.close();
    
    
    ofstream myfile3;
    myfile3.open ("valeurs_normales_3.txt");
    myfile3 << "Troisième tirage tirage tirage" << endl;
    for (int i=0; i<N; i++) {
        myfile3 << simulationNormale(3,i) << endl;
    }
    
    myfile3.close();
*/
    
    //-----------------------------------------
    // Calcul du prix et de l'intervalle de confiance
    //-----------------------------------------
 
    // optimisation
    vector<double> theta_3d_opt = newton_D3_panier(S0, r, sigma, T, K, N, simulationNormale, rho1, rho2, rho3,lambda);
    cout<<theta_3d_opt[0]<<" "<<theta_3d_opt[1]<<" "<<theta_3d_opt[2]<<endl;

    
    double prix=exp(-r*T)*mean_f_panier(S0,r,sigma,T,K,rho1,rho2,rho3,lambda,simulationNormale,theta_3d_opt);
    cout<<"Prix= "<<prix<< endl;
    
    double prixSansTetha=exp(-r*T)*mean_g_panier(S0,r,sigma,T,K,rho1,rho2,rho3, lambda,simulationNormale);
    cout<<"Prix sans theta = "<< prixSansTetha << endl;
    
    double intConfiance_panier=ic_f_panier(S0, r, sigma, T, K, rho1, rho2, rho3,lambda,simulationNormale,theta_3d_opt);
    cout <<"racine de la variance sur N= " << intConfiance_panier << endl;
    
    vector<double> intConfiance_X_Y=ic_f_panier_put(S0, r, sigma, T, K, rho1, rho2, rho3,lambda,simulationNormale);
    cout <<"variance de X-Y = " << intConfiance_X_Y[0] << " "<< "Variance de X " << intConfiance_X_Y[1] << endl;

    
    temps=clock();
    cout<<"temps nécessaire pour le calcul du prix et de son IC= "<< temps/(double)CLOCKS_PER_SEC << endl;
    
    
    ///ALTIPLANO

    /*
    // optimisation
    vector<double> theta_3d_opt = newton_D3_altiplano(S0, r, sigma, T, K, N, simulationNormale, rho1, rho2, rho3);
    cout<<theta_3d_opt[0]<<" "<<theta_3d_opt[1]<<" "<<theta_3d_opt[2]<<endl;
    
    double prix2=exp(-r*T)*mean_f_altiplano(S0,r,sigma,T,K,rho1,rho2,rho3,simulationNormale,theta_3d_opt);
    cout<<"Prix theta* = "<<prix2<< endl;
    
    double prix=exp(-r*T)*mean_g_altiplano(S0,r,sigma,T,K,rho1,rho2,rho3,simulationNormale);
    cout<<"Prix Altiplano= "<<prix<< endl;
    double intConfiance_altiplano=ic_f_altiplano(S0, r, sigma, T, K, rho1, rho2, rho3, simulationNormale,theta_3d_opt);
    cout <<"racine de la variance sur N (optimisé) = " << intConfiance_altiplano << endl;
    
    double intConfiance_altiplano_nul=ic_f_altiplano(S0, r, sigma, T, K, rho1, rho2, rho3, simulationNormale,theta_3d);
    cout <<"racine de la variance sur N (non optimisé) = " << intConfiance_altiplano_nul << endl;

    temps=clock();
    cout<<"temps nécessaire pour le calcul du prix et de son IC= "<< temps/(double)CLOCKS_PER_SEC << endl;
     */
    ///////---------------
    //
    //  Dernière question
    //
    /////////////---------
    
    // le premier terme c'est le prix d'un call sur panier, le second terme c'est le prix d'un put sur panier
    //double derniereQuestion = exp(-r*T)*(mean_g_panier(S0, r, sigma, T, K, rho1, rho2, rho3, lambda, simulationNormale)-mean_g_panier_put(S0, r, sigma, T, K, rho1, rho2, rho3, lambda, simulationNormale));
    //cout << derniereQuestion << endl;
    
    //-----------------------------------------------------------------------
    // Tests matrice
    //-----------------------------------------------------------------------
    /*
    matrice xxx(3, 3,0);
    xxx(3,3)=2;xxx(2,2)=1;xxx(1,1)=1;xxx(1,2)=3;xxx(2,1)=3;
    xxx+=xxx;
    vector<double> vec (3,1);
    vec[0]=0;
    vec[1]=0;
    vec[2]=4;
    
    
    
    int it=0;
    cout << xxx << endl;
    cout << xxx.inv_3d();
    vector<double> res=xxx.inv_3d()*vec;
    cout << res[0] << " " << res[1] << " " << res[2] << endl;
    */
    
    /*
    matrice R(3,3,0);
    R(1,1)=1;R(2,2)=1;R(3,3)=1;
    R(1,2)=0.5;R(1,3)=0.5;R(2,3)=0.5;
    R(2,1)=0.5;R(3,1)=0.5;R(3,2)=0.5;
    matrice chol=R.cholesky_3d();
    cout << chol;
    vector<double> premCol={chol(2,1),chol(2,2),chol(2,3)};
    cout << (chol*premCol)[0] << endl;
    cout << (chol*premCol)[1] << endl;
    cout << (chol*premCol)[2] << endl;
     */
    
    return 0;
}
