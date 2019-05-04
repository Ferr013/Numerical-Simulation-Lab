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
#include <sstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;
ofstream Eq;//Used to check

int main()
{ 
  Input(); //Inizialization
  Eq.open("out_eq.txt");
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  Eq.close();
  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");

   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> restart;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  if(restart==false){
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }else{
    ReadInput.open("config.final");
    cout << "Initial configuration set from config.final" << endl;
    for (int i=0; i<nspin; ++i)
    {
      ReadInput >> s[i];
    }
    ReadInput.close();
  }
    accepted=0;attempted=0;
    
  //Evaluate energy etc. of the initial configuration
    Measure();

  //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
      double acc = min(1.,exp(-beta*(Boltzmann(-s[o],o)-Boltzmann(s[o],o))));
      if(acc == 1 || acc>rnd.Rannyu()){
        s[o] = -s[o];
        accepted++;
      }
      attempted++;
    }
    else //Gibbs sampling
    {
      double probGibbs = 1./(1+exp(beta*(Boltzmann(1,o)-Boltzmann(-1,o))));
      if(probGibbs>rnd.Rannyu()){
        s[o] = +1;
      }else{
        s[o] = -1;
      }
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
  //Eq << u << endl; Used for thermal equilibrazation of the algorithm
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
    
    cout << "Block number " << iblk << endl;
    if(metro==1)cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    if(metro==1){
      stringstream ss;
      ss << "output.ene." << (int)((temp-0.5)/0.15);
      Ene.open ( ss.str().c_str(),ios::app);
      ss.str("");ss << "output.heat." << (int)((temp-0.5)/0.15);
      Heat.open(ss.str().c_str(),ios::app);
      ss.str("");ss << "output.mag." << (int)((temp-0.5)/0.15);
      Mag.open(ss.str().c_str(),ios::app);
      ss.str("");ss << "output.chi." << (int)((temp-0.5)/0.15);
      Chi.open(ss.str().c_str(),ios::app);
    }else{
      stringstream ss;
      ss << "outputG.ene." << (int)((temp-0.5)/0.15);
      Ene.open ( ss.str().c_str(),ios::app);
      ss.str("");ss << "outputG.heat." << (int)((temp-0.5)/0.15);
      Heat.open(ss.str().c_str(),ios::app);
      ss.str("");ss << "outputG.mag." << (int)((temp-0.5)/0.15);
      Mag.open(ss.str().c_str(),ios::app);
      ss.str("");ss << "outputG.chi." << (int)((temp-0.5)/0.15);
      Chi.open(ss.str().c_str(),ios::app);
    }
    if(h==0){
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene <<  iblk <<  " " << stima_u <<  " " << glob_av[iu]/(double)iblk <<  " " << err_u << endl;

    stima_c = beta*beta *((blk_av[ic]/blk_norm)-pow(stima_u*nspin,2.)); //Heat Capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat <<  iblk <<  " " << stima_c <<  " " << glob_av[ic]/(double)iblk/(double)nspin <<  " " << err_c/(double)nspin << endl;

    stima_x = beta*((blk_av[ix]/blk_norm/(double)nspin)-stima_m*stima_m); //Magnetic Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi <<  iblk <<  " " << stima_x <<  " " << glob_av[ix]/(double)iblk <<  " " << err_x << endl;
    }else{
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetisation
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag <<  iblk <<  " " << stima_m <<  " " << glob_av[im]/(double)iblk <<  " " << err_m << endl;
    }
    Ene.close();
    Heat.close();
    Mag.close();
    Chi.close();

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
