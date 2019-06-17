/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Giovanni Ferrami
_/    _/  _/_/_/  _/_/_/_/ email: g.ferrami@gmail.com
*****************************************************************
*****************************************************************/

#ifndef __GA_
#define __GA_

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;


//Random numbers
#include "random.h"
Random rnd;

double inline Random(){
    return rnd.Rannyu();
}

int randseed, wd=18;
uint gene_length=30;
double temp_init = 3;
int temperature_steps = 1000;
int steps_per_T = 10000;     
double temp;
int mode = 0; double A=0;
int accepted_1 = 0,accepted_2 = 0,accepted_3 = 0;


class City{
public:
    City();
    City(double x, double y, int i)
    {
        m_x = x;m_y = y;index = i;
    }
    ~City(){};

    double GetX(){return m_x;}
    double GetY(){return m_y;}
    int GetIndex(){return index;}

    void SetX(double x){m_x = x;}
    void SetY(double y){m_y = y;}
    void SetIndex(int i){index = i;}
private:
    double m_x,m_y;
    int index;
};

class Gene{
public:
    Gene() 
    { 
        cout << "--- Empty gene ---" << endl;  
    } 
    Gene(int l, vector<City>* c, bool shuffle) 
    { 
        length=l;
        for(uint i=0;i<length;i++){
            gene.push_back(c->at(i));
        }
        if(!CheckLength()) cout << "ERROR: gene length isn't what we want" << endl;
        if(shuffle) Shuffle();
        cost = Cost();
    } 
    ~Gene(){ 
        gene.clear(); 
    } 
    void AddElement(City city){
        gene.push_back(city);
    }
    City GetElement(int i){
        return gene[i];
    }
    uint GetLength(){
        length = gene.size();
        return length;
        }
    double GetCost(){
        return Cost();
        }
    bool CheckLength(){return (gene.size()==length);}

    void SwapCity(int x, int y){
        City temp = City(0,0,0);
        temp = gene[x];
        gene[x] = gene[y];
        gene[y] = temp;
    }

    double Cost(){
        double r = 0.;
        for (uint i = 0; i < GetLength(); i++)
        {   
            double j = i+1;
            if(i==length-1)j=0; //PBC
            r += sqrt(pow((gene[i].GetX()-gene[j].GetX()),2)+pow((gene[i].GetY()-gene[j].GetY()),2));
        }
        return r;
    }
    void Shuffle(){
        for (uint i=0; i<GetLength(); i++){
            //int index = (int)length*rnd.Rannyu();
            random_shuffle(gene.begin(), gene.end());
        }
    }
    void Print(){
        cout << "--------cities--------"<<endl;
        for (uint i=0; i<GetLength(); i++){
            cout << "City #" << gene[i].GetIndex() << " at "<<  gene[i].GetX() << " , " << gene[i].GetY() << endl;
        }
        cout << "Cost: " << Cost() << " $"<< endl;
    }
    bool CheckIsSalesman(){
        int check = 0;
        int theory = 0;
        for (uint i=0; i<GetLength(); i++){
            theory += pow(i,3);
            check += pow(gene[i].GetIndex(),3);
        }
        return check == theory;
    }

    void Mutation_Pair(){
        uint index = (uint)(GetLength()*Random());
        if(index != GetLength()-1) SwapCity(index,index+1);
        else SwapCity(index,0);
    }
    void Mutation_Shift(){
        int l = GetLength();
        for(int i=0;i<l-1;i++){
            SwapCity(i, i+1);
        }
    }
    void Mutation_Shift_N(){
        int m = 3;//METTERE IN FILE INPUT
        uint index = (uint)((GetLength()-m-1)*Random());
        for(int i=m;i>0;i--){
            SwapCity(index+i-1, index+i);
            SwapCity(index+i, index+i+1);
        }
    }
    void Mutation_Inversion(){
        int m = 3;//METTERE IN FILE INPUT
        uint index = (uint)((GetLength()-m-1)*Random()); 
        for(int i=0;i<=(int)m/2;i++){
            SwapCity(index+i, index+m-i);
        }
    }
    void Mutation_SwitchEnd(){
        int m = 5;//METTERE IN FILE INPUT 
        for(int i=0;i<=m;i++){
            SwapCity(i, GetLength()-m-1+i);
        }
    }
private:  
    uint length;
    vector<City> gene;
    double cost;
};



class Population{
public:
    Population();
    Population(uint pop_dim, uint gene_length){
        dim=pop_dim;
        for(uint i=0;i<dim;i++){
            Gene gene = Gene(gene_length, 0, true);
            population.push_back(gene);
        }
        //Print();
        Sort();
    }
    Population(uint pop_dim, uint gene_length, vector<City>* c){
        dim=pop_dim;
        for(uint i=0;i<dim;i++){
            Gene gene = Gene(gene_length, c, true);
            population.push_back(gene);
        }
        //Print();
        Sort();
    }
    ~Population(){ 
        population.clear(); 
    }
    void Clear(){ 
        population.clear(); 
    }
    void AddGene(Gene gene){
        population.push_back(gene);
    }
    Gene GetGene(int i){
        return population[i];
    } 
    uint GetDim(){
        dim = population.size();
        return dim;
    }
    static inline bool sort_by_cost (Gene& x, Gene& y)
    {
        return (x.GetCost() < y.GetCost());
    }
    void Sort(){
        sort (population.begin(), population.end(), sort_by_cost); 
    }

    /*
    void Mutate(double mutation_rate){
        for(uint i=0; i<GetDim(); i++){
            if(Random()<mutation_rate)Mutation_Pair(&population[i]);
            if(Random()<mutation_rate)Mutation_Shift(&population[i]);
            if(Random()<mutation_rate)Mutation_Shift_N(&population[i]);
            if(Random()<mutation_rate)Mutation_Inversion(&population[i]);
        }
    }
    */
    Gene SelectGene(double exp){//Select gene with higher probablity to lowest cost ones -> exp>1
        double s = Random(); //METTERE exp IN FILE INPUT
        double prob = pow(s,exp);
        double pos = (int)(dim*prob);
        return population[pos];

    }
    void Breed(Gene x, Gene y, Population* new_pop){
        uint length = x.GetLength();
        uint pos = (int)(length*Random());
        Gene a = Gene(0,0,false);Gene b = Gene(0,0,false);
        for(uint i=0;i<pos;i++){
            a.AddElement(x.GetElement(i));
            b.AddElement(y.GetElement(i));
        }
        while(a.GetLength()<length){
            for(uint i=0;i<length;i++){
                bool ok = true;
                for (uint j=0; j<a.GetLength(); j++)
                {
                    if(y.GetElement(i).GetIndex()==a.GetElement(j).GetIndex())ok=false;
                }
                if(ok){
                    a.AddElement(y.GetElement(i));
                }
            }
        }
        while(b.GetLength()<length){
            for(uint i=0;i<length;i++){
                bool ok = true;
                for (uint j=0; j<b.GetLength(); j++)
                {
                    if(x.GetElement(i).GetIndex()==b.GetElement(j).GetIndex())ok=false;
                }
                if(ok){
                    b.AddElement(x.GetElement(i));
                }
            }
        }
        new_pop->AddGene(a);
        new_pop->AddGene(b);
    }
    void SaveBest(int saved_per_generation, double selection_exp,Population* new_pop){
        for (int i = 0; i < saved_per_generation; i++)
        {
            Gene a = SelectGene(selection_exp);
            new_pop->AddGene(a);
        }
    }

    void Print(){
        cout << "--------Population--------"<<endl;
        for (uint i=0; i<GetDim(); i++){
            cout << "Gene #" << i << " | Cost: "  <<population[i].GetCost() << endl;
        }
    }
private:  
    uint dim;
    vector<Gene> population;
};



void Input(int rank, vector<City> *cities)
{
    ifstream ReadInput;

    if(rank==0)cout << "**************** SIMULATED ANNEALING ALGORITHM ******************" << endl;
    if(rank==0)cout << "           SOLVING THE TRAVELLING SALESMAN PROBLEM           " << endl;
    if(rank==0)cout << "             CODE ADAPTED TO PARALLEL COMPUTATION            " << endl << endl;

    //Read input informations
    
    ReadInput.open("input.dat");

    ReadInput >> randseed;
    //srand(randseed); //Used std random generator for the serial code
    if(rank==0)cout << "Random seed = " << randseed << endl;

    ReadInput >> gene_length;
    if(rank==0)cout << "Travel length = " << gene_length << endl;
        
    ReadInput >> temp_init;
    ReadInput >> temperature_steps;
    ReadInput >> steps_per_T;
    ReadInput >> mode;

    if(rank==0){
    cout << "Initial temperature = " << temp_init << endl;
    cout << "Temperature steps = " << temperature_steps << endl;
    cout << "Steps at fixed temperature = " << steps_per_T << endl;
    cout << "Mode= " << mode << endl << endl;}
    
    ReadInput.close();

    ////////////////////////////////////////////////Set the first Random Generator to create in every node the same city
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
    rnd.SaveSeed();

    ////////////////////////////////////////////////Create cities map
    if(mode == 0){   //Cities on a circumference
    double angle=2*M_PI/(double)gene_length;
    if(rank==0)cout << "Cities generated on a circumference "<< endl << endl;
    for (uint i=0; i<gene_length; i++){
        City city = City(cos(angle*i), sin(angle*i) ,i);
        cities->push_back(city);
    }}
    if(mode == 1){   //Cities inside a square
    double side = 1;
    if(rank==0)cout << "Cities generated inside a square "<< endl << endl;
    for (uint i=0; i<gene_length; i++){
        City city = City(Random()*side, Random()*side ,i);
        cities->push_back(city);
    }}

    /////////////////////////////////////////////////Initialize again the Random Generator to avoid that each node performs the same steps
    Primes.open("Primes");
    if (Primes.is_open()){
        //Read different primes for each node
        for(int i=0;i<rank+1;i++){
            Primes >> p1 >> p2 ;
        }
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    input.open("seed.in");
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
    rnd.SaveSeed();
}

#endif
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Giovanni Ferrami
_/    _/  _/_/_/  _/_/_/_/ email: g.ferrami@gmail.com
*****************************************************************
*****************************************************************/