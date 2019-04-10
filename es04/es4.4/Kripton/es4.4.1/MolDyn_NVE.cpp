/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int start_from_previous = 0;
int counter = 0;
const int iblock = 30;
const double arg_rho = 0.364, arg_T0=164,arg_m=83.798;
const double kb = 8.61673324*0.00001;

int main(){ 
  Input();             //Inizialization
  
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
        if(istep%(m_block*10)==0){
          //cout << "temp measured" << sum_temp/(double)m_block << endl;
          //cout << "Block finished at " << istep<< endl;
          block_epot[counter]=sum_pot/(double)m_block;
          block_ekin[counter]=sum_kin/(double)m_block;
          block_etot[counter]=sum_etot/(double)m_block;
          block_temp[counter]=sum_temp/(double)m_block;
          block_press[counter]=sum_press/(double)m_block;
          counter++;
          sum_pot=0;sum_kin=0;sum_etot=0;sum_temp=0;sum_press=0;
        }
     }
  }
  ofstream BPot("ave_epot.dat");	
	ofstream BKin("ave_ekin.dat");
	ofstream BTemp("ave_temp.dat");
	ofstream BTot("ave_etot.dat");
  ofstream BPress("ave_press.dat");
  for(int i=0;i<iblock;i++){
	    b2_epot[i]=block_epot[i]*block_epot[i];
	    b2_ekin[i]=block_ekin[i]*block_ekin[i];
	    b2_etot[i]=block_etot[i]*block_etot[i];
	    b2_temp[i]=block_temp[i]*block_temp[i];
      b2_press[i]=block_press[i]*block_press[i];
      for(int j=0;j<(i+1);j++){
        ave_epot[i]+=block_epot[j];
        av2_epot[i]+=b2_epot[j];
        ave_ekin[i]+=block_ekin[j];
        av2_ekin[i]+=b2_ekin[j];
        ave_etot[i]+=block_etot[j];
        av2_etot[i]+=b2_etot[j];
        ave_temp[i]+=block_temp[j];
        av2_temp[i]+=b2_temp[j];
        ave_press[i]+=block_press[j];
        av2_press[i]+=b2_press[j];
      }
      ave_epot[i]= ave_epot[i]/double(i+1);
      av2_epot[i]= av2_epot[i]/double(i+1);
      ave_ekin[i]= ave_ekin[i]/double(i+1);
      av2_ekin[i]= av2_ekin[i]/double(i+1);
      ave_temp[i]= ave_temp[i]/double(i+1);
      av2_temp[i]= av2_temp[i]/double(i+1);
      ave_etot[i]= ave_etot[i]/double(i+1);
      av2_etot[i]= av2_etot[i]/double(i+1);
      ave_press[i]= ave_press[i]/double(i+1);
      av2_press[i]= av2_press[i]/double(i+1);
      err_epot[i]=Error(ave_epot[i], av2_epot[i],i);
      err_ekin[i]=Error(ave_ekin[i], av2_ekin[i],i);
      err_etot[i]=Error(ave_etot[i], av2_etot[i],i);
      err_temp[i]=Error(ave_temp[i], av2_temp[i],i);
      err_press[i]=Error(ave_press[i], av2_press[i],i);
      BPot <<(i+1)<<" "<< ave_epot[i]/(arg_T0*kb) <<" "<< err_epot[i]/(arg_T0*kb) << endl;
      BKin <<(i+1)<<" "<< ave_ekin[i]/(arg_T0*kb)  <<" "<< err_ekin[i]/(arg_T0*kb) << endl;
      BTot <<(i+1)<<" "<< ave_etot[i]/(arg_T0*kb)  <<" "<< err_etot[i]/(arg_T0*kb) << endl;
      BTemp <<(i+1)<<" "<< ave_temp[i]/(arg_T0)  <<" "<< err_temp[i]/(arg_T0) << endl;
      BPress <<(i+1)<<" "<< (ave_press[i]/(arg_T0*kb))*pow(arg_rho,3)  <<" "<< (err_press[i]/(arg_T0*kb))*pow(arg_rho,3) << endl;
    }
  BPot.close();BKin.close();BTot.close();BTemp.close();BPress.close();
  ConfFinal();         //Write final configuration to restart
  ConfPrevious();      //Write the second to last position reached by the simulation

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input
  ReadInput >> start_from_previous;
  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables


//Read second to last configuration of the previous run to renormalize velocities to the target temperature
  if(start_from_previous==1){
    cout << "Read initial configuration from file config.0 " << endl << endl;
      ReadConf.open("config.final");
      for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
      }
      ReadConf.close();

    cout << "Read configuration to renormalize velocities from file config.prev " << endl << endl;
    ReadConf.open("config.prev");
    for (int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConf.close();

  //Find x(t+dt) with a single step of Verlet algorithm and then evaluate the velocity
    double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
    for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }
    double sumv[3] = {0.0, 0.0, 0.0};
    for(int i=0; i<npart; ++i){ //Verlet integration scheme

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    //Renormalize velocity to match the target temperature
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = x[i] - vx[i] * delta;
      yold[i] = y[i] - vy[i] * delta;
      zold[i] = z[i] - vz[i] * delta;
    }
    return;

  }else{
    //Read initial configuration
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();
  //Prepare initial velocities
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
      vx[i] = (double)rand()/RAND_MAX - 0.5;
      vy[i] = (double)rand()/RAND_MAX - 0.5;
      vz[i] = (double)rand()/RAND_MAX - 0.5;

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = x[i] - vx[i] * delta;
      yold[i] = y[i] - vy[i] * delta;
      zold[i] = z[i] - vz[i] * delta;
    }
    return;
  }
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij, p, pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Press.open("output_press.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0; p = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       pij = (48.0/pow(dr,12) - 24.0/pow(dr,6));

//Potential energy and pressure
       v += vij; p += pij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy
    stima_press = rho*stima_temp + (1.0/(3*vol))*(v/(double)npart);

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_press << endl;

    sum_pot+=stima_pot;
    sum_kin+=stima_kin;
    sum_etot+=stima_etot;
    sum_temp+=stima_temp;
    sum_press+=stima_press;
    //cout << "temp " << stima_temp<< endl;
    //cout << "temp_sum " << sum_temp<< endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfPrevious(void){  //Write second to last position of the run
  ofstream WriteConf;

  cout << "Print second to last configuration to file config.prev " << endl << endl;
  WriteConf.open("config.prev");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;

}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double av, double av2, int n){
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
