/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Giovanni Ferrami
_/    _/  _/_/_/  _/_/_/_/ email: g.ferrami@gmail.com
*****************************************************************
*****************************************************************/
#include "GA_TSP.h"

using namespace std;

void Input(void);

int randseed, wd=18;
uint dim = 900, gene_length=30;//METTERE IN FILE INPUT
uint saved_per_generation = 15; //METTERE IN FILE INPUT
double selection_exp = 5; //METTERE IN FILE INPUT
double breeding_rate = 0.8; //METTERE IN FILE INPUT
double mutation_rate = 0.1; //METTERE IN FILE INPUT
int num_generations = 100; //METTERE IN FILE INPUT
int mode = 0;
vector<City> cities;

int main(){
    //cout << cities.size() << endl;
    Input(); //Inizialization
    //cout << cities.size() << endl;

    Population pop = Population(dim, gene_length, &cities);
    Population next_pop = Population(0, gene_length, &cities);

    //pop.GetGene(0).Print();

    ofstream output, best, first;
    if(mode == 0)output.open("evolution.dat");
    if(mode == 1)output.open("evolution_sq.dat");
    if(mode == 0)first.open("first.gene");
    if(mode == 1)first.open("first_sq.gene");
    Gene g = pop.GetGene(0);
    for (uint i = 0; i < g.GetLength(); i++){    
        City c = g.GetElement(i);
        first << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;
    }
    first.close();
    for (int i = 0; i < num_generations; i++)
    {
        //cout << "------- Generation #" << i << " -------" <<endl;   
        //cout << "Saving Best" << endl;
        pop.SaveBest(saved_per_generation,selection_exp,&next_pop);
        //cout << "Mutation" << endl;
        pop.Mutate(mutation_rate);
        //pop.Print();
        pop.Sort();
        //next_pop.Print();
        //cout << "Breeding" << endl;
        if((dim-saved_per_generation)%2!=0){
            Gene x = pop.SelectGene(selection_exp);
            next_pop.AddGene(x);
        }
        for(uint j = 0; j < (uint)(dim-saved_per_generation)/2; j++){
            Gene x = pop.SelectGene(selection_exp);
            Gene y = pop.SelectGene(selection_exp);
            if(Random()<breeding_rate){pop.Breed(x,y, &next_pop);}
            else{next_pop.AddGene(x);next_pop.AddGene(y);}
        }
        //cout << "After Breeding" << endl;
        next_pop.Sort();
        pop.Clear();
        for (uint j = 0; j < next_pop.GetDim(); j++){
            pop.AddGene(next_pop.GetGene(j));
        }
        //pop.Print();
        pop.Mutate(mutation_rate);
        pop.Sort();
        //pop.Data() write data on a log file
        cout << "Generation #" << i << " best Salesman Cost: " << pop.GetGene(0).GetCost() <<endl;
        output << setw(wd) << i <<  setw(wd) << pop.GetGene(0).GetCost() << endl;
        next_pop.Clear();
        //pop.Print();
    }
    output.close();
    if(mode == 0)best.open("best.gene");
    if(mode == 1)best.open("best_sq.gene");
    g = pop.GetGene(0);
    for (uint i = 0; i < g.GetLength(); i++){    
        City c = g.GetElement(i);
        best << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;
    }
    best.close();
    return 0;
}

void Input(void)
{
    ifstream ReadInput;

    cout << "**************** GENETIC ALGORITHM ******************" << endl;
    cout << "      SOLVING THE TRAVELLING SALESMAN PROBLEM        " << endl << endl;

    //Read input informations
    ReadInput.open("input.dat");

    ReadInput >> randseed;
    srand(randseed);
    cout << "Random seed = " << randseed << endl;

    ReadInput >> dim;
    cout << "Population dimensions = " << dim << endl;

    ReadInput >> gene_length;
    cout << "Gene length = " << gene_length << endl;

    ReadInput >> num_generations;
    cout << "Number generations = " << num_generations << endl << endl;
        
    ReadInput >> saved_per_generation;
    ReadInput >> selection_exp;
    ReadInput >> breeding_rate;
    ReadInput >> mutation_rate;
    ReadInput >> mode;

    cout << "Selection exponent = " << selection_exp << endl;
    cout << "Mutation rate = " << mutation_rate << endl;
    cout << "Breeding rate = " << breeding_rate << endl;
    cout << "Elite saved per generation = " << saved_per_generation << endl << endl;
    ReadInput.close();

    ////////////////////////////////////////////////Create cities map
    if(mode == 0){   //Cities on a circumference
    double angle=2*M_PI/(double)gene_length;
    cout << "Cities generated on a circumference "<< endl << endl;
    for (uint i=0; i<gene_length; i++){
        City city = City(cos(angle*i), sin(angle*i) ,i);
        cities.push_back(city);
    }}
    if(mode == 1){   //Cities inside a square
    double side = 1;
    cout << "Cities generated inside a square "<< endl << endl;
    for (uint i=0; i<gene_length; i++){
        City city = City(Random()*side, Random()*side ,i);
        cities.push_back(city);
    }}

}


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Giovanni Ferrami
_/    _/  _/_/_/  _/_/_/_/ email: g.ferrami@gmail.com
*****************************************************************
*****************************************************************/