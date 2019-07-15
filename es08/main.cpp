/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
double error(double, double, int);
double Psi(double, double, double);
double Psi_2(double, double, double);
double V(double);

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


//es5.1 from here ---------------------------------------------------------
    int M=100000; int N=100; int L=M/N;
	double ave[N],av2[N],sum_prog[N],su2_prog[N],err_prog[N];
    double step = 2.2;
    double mu = 0.78;
    double sigma = 0.6;
    int searching_parameters = 1; //DA USARE PER STIMARE MU E SIGMA, METTERE A ZERO PER AVERE DATA BLOCKING SUI VALORI MIGLIORI

    int nbins = 100;
    int istogrammi[nbins][N];   //questa matrice contiene in ogni riga l'istogramma riferito ad un singolo blocco
    
    for(int i=0; i<nbins; i++){
        for(int j=0; j<N; j++){
            istogrammi[i][j] = 0;
        }
    }

    for(int i=0; i<N;i++){
		ave[i]=0;av2[i]=0;sum_prog[i]=0;su2_prog[i]=0;err_prog[i]=0;
	}
    if (searching_parameters == 1){
        cout << "-------- Optimizing parameters mu and sigma --------" << endl;
        double min_energy = 99999.9;

        for (double mu_test = 0.77; mu_test <= 0.79; mu_test= mu_test+0.005){
            for (double sigma_test = 0.58; sigma_test <= 0.62; sigma_test= sigma_test+0.005){
                double mean = 0; double x = 0.; double x_new = 0.;
                int attempts=0, accepted=0;
                for(int i=0; i<M;i++){
                    x_new = x + (rnd.Rannyu()-0.5)*2*step;
                    double ratio = pow(Psi(x_new, mu_test, sigma_test),2)/pow(Psi(x, mu_test, sigma_test),2);
                    double acc = min(1.,ratio);
                    if(acc == 1 || acc>rnd.Rannyu()){
                        x = x_new;
                        accepted++;
                    }
                    attempts++;
                    mean += (-Psi_2(x,mu,sigma)*0.5 + V(x)*Psi(x,mu,sigma))/(Psi(x,mu,sigma));
                }
                mean = mean/(double)M;
                cout << "Energy : " << mean << endl;
                if(mean<min_energy){
                    min_energy=mean;
                    mu = mu_test; sigma = sigma_test;
                }
                cout << "Acceptance : " << (double)(accepted)/(double)(attempts) << " %"<< endl;
            }
        }
        cout << "Best values for " << endl;
        cout << "   Energy: " << min_energy << endl;
        cout << "   Mu: " << mu << endl;
        cout << "   Sigma: " << sigma << endl;
    }
    cout << "Evaluating energy with data blocking " << endl;
    ofstream out1("Energy.dat");
    ofstream out("Hist.dat");
    int attempts=0, accepted=0;
    double x,x_new,mean;
    for(int j = 0; j < N; j++)
    {
        mean = 0;
        x = 0.;
        for(int i=0; i<L;i++){
            x_new = x + (rnd.Rannyu()-0.5)*2*step;
            double ratio = pow(Psi(x_new, mu, sigma),2)/pow(Psi(x, mu, sigma),2);
            double acc = min(1.,ratio);
            if(acc == 1 || acc>rnd.Rannyu()){
                x = x_new;
                accepted++;
            }
            attempts++;
            mean += (-Psi_2(x,mu,sigma)*0.5 + V(x)*Psi(x,mu,sigma))/(Psi(x,mu,sigma));

            int index;
            index = int((x+2.5)/0.05);
            istogrammi[index][j]++;
        }
        cout << "Acceptance : " << (double)(accepted)/(double)(attempts) << " %"<< endl;

        ave[j] = mean/double(L);
        av2[j] = ave[j]*ave[j];
    }
    
    for(int i=0;i<N;i++){
        for(int j=0;j<(i+1);j++){
            sum_prog[i]+=ave[j];
            su2_prog[i]+=av2[j];
        }
        sum_prog[i]= sum_prog[i]/double(i+1);
        su2_prog[i]= su2_prog[i]/double(i+1);
        err_prog[i]=error(sum_prog[i], su2_prog[i],i);

        out1 <<(i+1)<<" "<< sum_prog[i] <<" "<< err_prog[i] << endl;
    }
    cout << "Energy : " <<  sum_prog[N-1]<< endl;
    out1.close();
    out1.clear();
    
    for(int i=0; i<N; i++){
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }

    cout << "HISTOGRAM DATA BLOCKING" << endl << endl;
    
    for(int ibin=0; ibin<nbins; ibin++){
        for(int i=0; i<N; i++){
            sum_prog[i] = 0;
            su2_prog[i] = 0;
            err_prog[i] = 0;
        }
        for(int i=0; i<N; i++){
            for(int j=0; j<i+1; j++){
                sum_prog[i] = sum_prog[i] + istogrammi[ibin][j];
                su2_prog[i] = su2_prog[i] + istogrammi[ibin][j]*istogrammi[ibin][j];
            }
            
            sum_prog[i] = sum_prog[i]/(double)(i+1);
            su2_prog[i] = su2_prog[i]/(double)(i+1);
            err_prog[i] = error(sum_prog[i],su2_prog[i],i);
            if(i == N-1){
                out << (-2.5+0.025)+(ibin*0.05) <<"   "<< sum_prog[i]/(double(M/N)*0.05) <<"   "<< err_prog[i]/(double(M/N)*0.05) << endl;
            }
        }
    }
    out.close();
    out.clear();
    rnd.SaveSeed();
    return 0;
}


double error(double av, double av2, int n){
	double err =0;
	if(n==0)return 0.;
	else 
	{
		err=sqrt(1./n*(av2-av*av));
	}
	return err;
}

double Psi(double x, double mu, double sigma){
    double den = 2*sigma*sigma;
    double r = exp(-pow(x-mu,2)/den) + exp(-pow(x+mu,2)/den);
    return r;
}
double Psi_2(double x, double mu, double sigma){
    double den = 2*sigma*sigma;
    double b = 1./pow(sigma,2);
    double a1 = pow((x-mu)/(sigma*sigma),2);
    double a2 = pow((x+mu)/(sigma*sigma),2);
    double r = exp(-pow(x-mu,2)/den)*(a1-b) + exp(-pow(x+mu,2)/den)*(a2-b) ;
    return r;
}
double V(double x){
    double r = pow(x,4)-2.5*pow(x,2);
    return r;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/