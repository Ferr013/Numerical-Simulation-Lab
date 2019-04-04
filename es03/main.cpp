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
#include <vector>
#include "random.h"

double GMB(double ,double , double , double , double, double);
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
//es3.1 from here ---------------------------------------------------------
	int M=100000; int N=100; int L=M/N;
	double S_0 = 100.; double T =1.; double K=100.; double r=0.1;
	double sigma =0.25; double S = 0.;
	vector<double> ave;vector<double> av2;vector<double> sum_prog;
	vector<double> su2_prog;vector<double> err_prog;
	ofstream out1("direct_call.dat");	
	ofstream out2("discrete_call.dat");
	ofstream out3("direct_put.dat");
	ofstream out4("discrete_put.dat");

//sampling directly the final asset price
   double C = 0;double P=0;
   for(int i = 0; i < N; i++){
      for(int j = 0; j < L; j++)
      {
         S=GMB(S_0,r,sigma,0,T,rnd.Gauss(0.,1.));
         C += max(0., S-K);
      }
      ave.push_back(exp(-r*T)*C/double(L));
      av2.push_back(ave[i]*ave[i]);
      C = 0;
   }
   for(int i=0;i<N;i++){
      double a1=0;double a2=0;
		for(int j=0;j<(i+1);j++){
			a1+=ave[j];a2+=av2[j];
		}
      sum_prog.push_back(a1/double(i+1));
		su2_prog.push_back(a2/double(i+1));
		err_prog.push_back(error(sum_prog[i], su2_prog[i],i));

		out1 <<(i+1)<<" "<< sum_prog[i] <<" "<< err_prog[i] << endl;
	}
   ave.clear();av2.clear();sum_prog.clear();su2_prog.clear();err_prog.clear();
   
   //riscrivi i vettori e risparmia cicli
   for(int i = 0; i < N; i++){
      for(int j = 0; j < L; j++)
      {
         S=GMB(S_0,r,sigma,0,T,rnd.Gauss(0.,1.));
         P += max(0., K-S);
      }
      ave.push_back(exp(-r*T)*P/double(L));
      av2.push_back(ave[i]*ave[i]);  P=0;
   }
   for(int i=0;i<N;i++){
		double a1=0;double a2=0;
		for(int j=0;j<(i+1);j++){
			a1+=ave[j];a2+=av2[j];
		}
      sum_prog.push_back(a1/double(i+1));
		su2_prog.push_back(a2/double(i+1));
		err_prog.push_back(error(sum_prog[i], su2_prog[i],i));

		out3 <<(i+1)<<" "<< sum_prog[i] <<" "<< err_prog[i] << endl;
	}
	out1.close();
	out3.close();
   ave.clear();av2.clear();sum_prog.clear();su2_prog.clear();err_prog.clear();


//sampling discretized GBM path of the asset price
   for(int i = 0; i < N; i++){
      for(int j = 0; j < L; j++)
      {
         double s=S_0; double t_0 = 0;double t_1 = 0;
         for(int k = 0; k < 100; k++)
         {  
            t_1 += 0.01;
            S=GMB(s,r,sigma,t_0,t_1,rnd.Gauss(0.,1.));
            t_0 = t_1; s = S;
         }
         C += max(0., S-K);
      }
      ave.push_back(exp(-r*T)*C/double(L));
      av2.push_back(ave[i]*ave[i]);
      C = 0;
   }
   for(int i=0;i<N;i++){
      double a1=0;double a2=0;
		for(int j=0;j<(i+1);j++){
			a1+=ave[j];a2+=av2[j];
		}
      sum_prog.push_back(a1/double(i+1));
		su2_prog.push_back(a2/double(i+1));
		err_prog.push_back(error(sum_prog[i], su2_prog[i],i));

		out2 <<(i+1)<<" "<< sum_prog[i] <<" "<< err_prog[i] << endl;
	}
   ave.clear();av2.clear();sum_prog.clear();su2_prog.clear();err_prog.clear();
   
   for(int i = 0; i < N; i++){
      for(int j = 0; j < L; j++)
      {
         double s=S_0; double t_0 = 0;double t_1 = 0;
         for(int k = 0; k < 100; k++)
         {  
            t_1 += 0.01;
            S=GMB(s,r,sigma,t_0,t_1,rnd.Gauss(0.,1.));
            t_0 = t_1; s = S;
         }
         P += max(0., K-S);
      }
      ave.push_back(exp(-r*T)*P/double(L));
      av2.push_back(ave[i]*ave[i]);
      P = 0;
   }
   for(int i=0;i<N;i++){
      double a1=0;double a2=0;
		for(int j=0;j<(i+1);j++){
			a1+=ave[j];a2+=av2[j];
		}
      sum_prog.push_back(a1/double(i+1));
		su2_prog.push_back(a2/double(i+1));
		err_prog.push_back(error(sum_prog[i], su2_prog[i],i));

		out4 <<(i+1)<<" "<< sum_prog[i] <<" "<< err_prog[i] << endl;
	}
	out2.close();
	out4.close();

   rnd.SaveSeed();
   return 0;
}

double GMB(double S_0,double mu, double sigma, double t0, double t1,double g){
	double S=S_0*exp((mu-sigma*sigma/2.)*(t1-t0)+sigma*g*sqrt(t1-t0));
   return S;
}

double error(double av, double av2, int n){
	double err =0;
	if(n==0)return 0.;
	else 
	{
		//cout << abs(av2-av*av) << endl;
		err=sqrt(1./double(n)*(av2-av*av));
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
