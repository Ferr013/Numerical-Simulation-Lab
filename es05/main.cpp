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
double S(double, double, double);
double P(double, double, double);

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

	ofstream out1("1s_ave.dat");
    ofstream out2("2p_ave.dat");

    for(int i=0; i<N;i++){
		ave[i]=0;av2[i]=0;sum_prog[i]=0;su2_prog[i]=0;err_prog[i]=0;
	}

    double x,y,z,x_new,y_new,z_new,mean;
    for(int j = 0; j < N; j++)
    {
        mean = 0;
        //condizione iniziale per 1s è l'origine
        double step = 1.2;
        x = 0.; y = 0. ; z = 0.;
        for(int i=0; i<L;i++){
            x_new = x + (rnd.Rannyu()-0.5)*2*step;
            y_new = y + (rnd.Rannyu()-0.5)*2*step;
            z_new = z + (rnd.Rannyu()-0.5)*2*step;

            double acc = min(1.,S(x_new,y_new,z_new)/S(x,y,z));
            if(acc == 1 || acc>rnd.Rannyu()){
                x = x_new; y = y_new; z = z_new;
            }
            mean += sqrt(x*x + y*y + z*z);
        }
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
	out1.close();
    
    for(int i=0; i<N;i++){
		ave[i]=0;av2[i]=0;sum_prog[i]=0;su2_prog[i]=0;err_prog[i]=0;
	}

    for(int j = 0; j < N; j++)
    {
        mean = 0;
        //condizione iniziale per 2p è
        double step = 3.;
        x = 0.; y = 0. ; z = 2.5;
        for(int i=0; i<L;i++){
            x_new = x + (rnd.Rannyu()-0.5)*2*step;
            y_new = y + (rnd.Rannyu()-0.5)*2*step;
            z_new = z + (rnd.Rannyu()-0.5)*2*step;

            double acc = min(1.,P(x_new,y_new,z_new)/P(x,y,z));
            if(acc == 1 || acc>rnd.Rannyu()){
                x = x_new; y = y_new; z = z_new;
            }
            mean += sqrt(x*x + y*y + z*z);
        }
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

		out2 <<(i+1)<<" "<< sum_prog[i] <<" "<< err_prog[i] << endl;
	}
	out2.close();
    
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

double S(double x, double y, double z){
    double r = sqrt(x*x + y*y + z*z);
    return pow(pow(M_PI,-0.5)*exp(-r),2);
}
double P(double x, double y, double z){
    double r = sqrt(x*x + y*y + z*z);
    return pow(pow(M_PI/2.,-0.5)/8.*r*exp(-r/2.)*z/r,2);
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