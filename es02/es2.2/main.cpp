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


//es2.2.1 from here ---------------------------------------------------------
   int N=10000; int L=100;
	double ave[N],av2[N],sum_prog[N],su2_prog[N],err_prog[N];
	ofstream out1("es2_2_1.dat");
    ofstream out11("es2_2_1_step.dat");	

	for(int i=0; i<N;i++){
		ave[i]=0;av2[i]=0;sum_prog[i]=0;su2_prog[i]=0;err_prog[i]=0;
	}

	for(int i=0; i<N;i++){
        int x_tot = 0; int y_tot = 0; int z_tot = 0;
		for(int j=0; j<L; j++){
            int r = rnd.Rannyu(0,6);
            if(r==0) x_tot++;
            if(r==1) x_tot--;
            if(r==2) y_tot++;
            if(r==3) y_tot--;
            if(r==4) z_tot++;
            if(r==5) z_tot--;

            if(i==0){
            out11 << j << " "<<sqrt(x_tot*x_tot + y_tot*y_tot + z_tot*z_tot) << endl;
            }
		}
		ave[i]=sqrt(x_tot*x_tot + y_tot*y_tot + z_tot*z_tot);
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

		out1 <<(i+1)<<" "<< sum_prog[i] <<" "<< err_prog[i] << endl;
	}
	out1.close();
    out11.close();


//The correct solution here
    ofstream out111("es2_2_1_sm.dat");
    double rw[L],rw2[L],sum_rw[L],su2_rw[L],err_rw[L];
    for(int i=0; i<L;i++){
		rw[i]=0;rw2[i]=0;sum_rw[i]=0;su2_rw[i]=0;err_rw[i]=0;
	} 
    for(int i=0; i<N;i++){
        int x_tot = 0; int y_tot = 0; int z_tot = 0;
		for(int j=0; j<L; j++){
            int r = rnd.Rannyu(0,6);
            if(r==0) x_tot++;
            if(r==1) x_tot--;
            if(r==2) y_tot++;
            if(r==3) y_tot--;
            if(r==4) z_tot++;
            if(r==5) z_tot--;
            rw[j]+= sqrt(x_tot*x_tot + y_tot*y_tot + z_tot*z_tot);
		}
    }	

    for(int i=0;i<L;i++){
        rw[i]=rw[i]/double(N);
	    rw2[i]=rw[i]*rw[i];
        for(int j=0;j<(i+1);j++){
			sum_rw[i]+=rw[j];
			su2_rw[i]+=rw2[j];
		}
		sum_rw[i]= sum_rw[i]/double(i+1);
		su2_rw[i]= su2_rw[i]/double(i+1);
		err_rw[i]=error(sum_rw[i], su2_rw[i],i);
        out111 <<(i+1)<<" "<< rw[i] <<" "<< err_rw[i] << endl;
    }
    out111.close();

//es2.2.2 from here -------------------------

    ofstream out2("es2_2_2.dat");	
    ofstream out22("es2_2_2_step.dat");

	for(int i=0; i<N;i++){
		ave[i]=0;av2[i]=0;sum_prog[i]=0;su2_prog[i]=0;err_prog[i]=0;
	}

	
	for(int i=0; i<N;i++){
        double x_tot = 0; double y_tot = 0; double z_tot = 0;
		for(int j=0; j<L; j++){
            double phi = rnd.Rannyu(0,2*M_PI);
            double theta = acos(1-2*rnd.Rannyu());
            
            x_tot += sin(theta)*cos(phi);
            y_tot += sin(theta)*sin(phi);
            z_tot += cos(theta);

            if(i==0){
            out22 << j << " "<<sqrt(x_tot*x_tot + y_tot*y_tot + z_tot*z_tot) << endl;
            }
		}
		ave[i]=sqrt(x_tot*x_tot + y_tot*y_tot + z_tot*z_tot);
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
    out22.close();



//The correct solution here
    ofstream out222("es2_2_2_sm.dat");
    for(int i=0; i<L;i++){
		rw[i]=0;rw2[i]=0;sum_rw[i]=0;su2_rw[i]=0;err_rw[i]=0;
	}
    for(int i=0; i<N;i++){
        double x_tot = 0; double y_tot = 0; double z_tot = 0;
		for(int j=0; j<L; j++){
            double phi = rnd.Rannyu(0,2*M_PI);
            double theta = acos(1-2*rnd.Rannyu());
            x_tot += sin(theta)*cos(phi);
            y_tot += sin(theta)*sin(phi);
            z_tot += cos(theta);

            rw[j]+= sqrt(x_tot*x_tot + y_tot*y_tot + z_tot*z_tot);
		}
    }	

    for(int i=0;i<L;i++){
        rw[i]=rw[i]/double(N);
	    rw2[i]=rw[i]*rw[i];
        for(int j=0;j<(i+1);j++){
			sum_rw[i]+=rw[j];
			su2_rw[i]+=rw2[j];
		}
		sum_rw[i]= sum_rw[i]/double(i+1);
		su2_rw[i]= su2_rw[i]/double(i+1);
		err_rw[i]=error(sum_rw[i], su2_rw[i],i);
        out222 <<(i+1)<<" "<< rw[i] <<" "<< err_rw[i] << endl;
    }
    out222.close();
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