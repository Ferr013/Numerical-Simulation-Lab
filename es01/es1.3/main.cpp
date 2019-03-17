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

double error(double, double, int);

using namespace std;
 
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
/*
   for(int i=0; i<20; i++){
      cout << rnd.Rannyu() << endl;
   }
*/	
//es1.3 from here ---------------------------------------------------------
	double l=1.; double d=1.3;
   int M=500000; int N=100; int L=M/N;
	double ave[N],av2[N],sum_prog[N],su2_prog[N],err_prog[N];
	ofstream out1("es1_3.dat");	

	for(int i=0; i<N;i++){
		ave[i]=0;av2[i]=0;sum_prog[i]=0;su2_prog[i]=0;err_prog[i]=0;
	}

	for(int i=0; i<N;i++){
		double hits = 0;
		for(int j=0; j<L; j++){
         double angle = atan(rnd.Rannyu()/rnd.Rannyu());
			if(d*rnd.Rannyu()-l*sin(angle)<0) hits+=1.;
		}
		ave[i]=2*l*double(L)/(hits*d);
		av2[i]=ave[i]*ave[i];	
	}

	for(int i=0;i<N;i++){
		for(int j=0;j<(i+1);j++){
			sum_prog[i]+=ave[j];
			su2_prog[i]+=av2[j];
		}
		sum_prog[i]= sum_prog[i]/double(i+1);
		su2_prog[i]= su2_prog[i]/double(i+1);
		err_prog[i]=error(sum_prog[i], su2_prog[i],i);

		out1 <<(i+1)*L<<" "<< sum_prog[i] <<" "<< err_prog[i] << endl;
	}
	out1.close();


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

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
