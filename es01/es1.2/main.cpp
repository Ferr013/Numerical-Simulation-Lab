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
//es1.2 from here ---------------------------------------------------------
	int M=10000;
   ofstream out0("std_N_1.dat");
	ofstream out1("exp_N_1.dat");
	ofstream out2("lor_N_1.dat");
	for(int i = 0; i < M; i++)
	{  
      int N=1;
      double sums=0.;double sume=0.;double suml=0.;
      for(int j = 0; j < N; j++)
      {  
         sums += rnd.Rannyu();
         sume += rnd.Exp(1.);
         suml += rnd.Lorentz(0.,1.);
      }
      out0<< i << ' ' << sums/N << endl;
		out1<< i << ' ' << sume/N << endl;
		out2<< i << ' ' << suml/N << endl;
	}
	out0.close();out1.close();out2.close();

   ofstream out02("std_N_2.dat");
	ofstream out12("exp_N_2.dat");
	ofstream out22("lor_N_2.dat");
	for(int i = 0; i < M; i++)
	{  
      int N=2;
      double sums=0.;double sume=0.;double suml=0.;
      for(int j = 0; j < N; j++)
      {  
         sums += rnd.Rannyu();
         sume += rnd.Exp(1.);
         suml += rnd.Lorentz(0.,1.);
      }
      out02<< i << ' ' << sums/N << endl;
		out12<< i << ' ' << sume/N << endl;
		out22<< i << ' ' << suml/N << endl;
	}
	out02.close();out12.close();out22.close();

   ofstream out00("std_N_10.dat");
	ofstream out10("exp_N_10.dat");
	ofstream out20("lor_N_10.dat");
	for(int i = 0; i < M; i++)
	{  
      int N=10;
      double sums=0.;double sume=0.;double suml=0.;
      for(int j = 0; j < N; j++)
      {  
         sums += rnd.Rannyu();
         sume += rnd.Exp(1.);
         suml += rnd.Lorentz(0.,1.);
      }
      out00<< i << ' ' << sums/N << endl;
		out10<< i << ' ' << sume/N << endl;
		out20<< i << ' ' << suml/N << endl;
	}
	out00.close();out10.close();out20.close();

   ofstream out000("std_N_100.dat");
	ofstream out100("exp_N_100.dat");
	ofstream out200("lor_N_100.dat");
	for(int i = 0; i < M; i++)
	{  
      int N=100;
      double sums=0.;double sume=0.;double suml=0.;
      for(int j = 0; j < N; j++)
      {  
         sums += rnd.Rannyu();
         sume += rnd.Exp(1.);
         suml += rnd.Lorentz(0.,1.);
      }
      out000<< i << ' ' << sums/N << endl;
		out100<< i << ' ' << sume/N << endl;
		out200<< i << ' ' << suml/N << endl;
	}
	out000.close();out100.close();out200.close(); 



   rnd.SaveSeed();
   return 0;
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
