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
//es1.1 from here ---------------------------------------------------------
	int M=100000; int N=100; int L=M/N;
	double ave[N],av2[N],sum_prog[N],su2_prog[N],err_prog[N];
	ofstream out1("es1_1_mean.dat");	

	for(int i=0; i<N;i++){
		ave[i]=0;av2[i]=0;sum_prog[i]=0;su2_prog[i]=0;err_prog[i]=0;
	}

	for(int i=0; i<N;i++){
		double sum = 0;
		for(int j=0; j<L; j++){
			sum+=rnd.Rannyu();
		}
		ave[i]=sum/L;
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
//es1.1.2--------------------------------------
	ofstream out2("es1_1_var.dat");	
	for(int i=0; i<N;i++){
		ave[i]=0;av2[i]=0;sum_prog[i]=0;su2_prog[i]=0;err_prog[i]=0;
	}

	for(int i=0; i<N;i++){
		double sum = 0;
		for(int j=0; j<L; j++){
			sum+=pow((rnd.Rannyu()-0.5),2);
		}
		ave[i]=sum/L;
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

		out2 <<(i+1)*L<<" "<< sum_prog[i] <<" "<< err_prog[i] << endl;
	}

	out2.close();
//es1.1.3-----------------------------------
	ofstream out3("es1_1_chi.dat");	
	M=100;int n=10000;
	int r=n/M;
	double d[100];
	double chi=0;
	for(int j=0; j<M; j++){
			d[j]=0;
		}

	for(int i=0;i<M;i++){
		int l=0;chi=0;
		for(int j=0;j<n;j++){
			l = int(100*rnd.Rannyu());
			d[l]+=1;
		}
		for(int j=0; j<M; j++){
			chi += (d[j]-r)*(d[j]-r)/r;
			d[j]=0;
		}
		out3<<i+1<<" "<<chi<<endl;		
	}
	out3.close();
   rnd.SaveSeed();
   return 0;
}

double error(double av, double av2, int n){
	double err =0;
	if(n==0)return 0.;
	else 
	{
		//cout << abs(av2-av*av) << endl;
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
