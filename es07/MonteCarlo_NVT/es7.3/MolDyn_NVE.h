/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie,igofr;
double bin_size,nbins,stima_gofr,stima_pot, stima_kin, stima_etot, stima_temp, stima_press;
double sum_pot=0, sum_kin=0, sum_etot=0, sum_temp=0, sum_press=0; 
double walker[1000];

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];
const int m_block=100;
double block_epot[m_block],block_ekin[m_block],block_etot[m_block],block_temp[m_block],block_press[m_block];
double b2_epot[m_block],b2_ekin[m_block],b2_etot[m_block],b2_temp[m_block],b2_press[m_block];
double ave_epot[m_block],ave_ekin[m_block],ave_etot[m_block],ave_temp[m_block],ave_press[m_block];
double av2_epot[m_block],av2_ekin[m_block],av2_etot[m_block],av2_temp[m_block],av2_press[m_block];
double err_epot[m_block],err_ekin[m_block],err_etot[m_block],err_temp[m_block],err_press[m_block];
double blk_norm,err_gdir;
double r, gdir;
double blk_av[1000];
double glob_av[1000],glob_av2[1000];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfPrevious(void);
void ConfXYZ(int);
void Measure(void);
void Accumulate(void);
void GetGave(void);
void Results(void);
double Force(int, int);
double Pbc(double);
double Error(double, double , int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
