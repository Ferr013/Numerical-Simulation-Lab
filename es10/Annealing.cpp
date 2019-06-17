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

vector<City> cities;

int main(){
    Input(&cities);
    ofstream output, best, first;

    Gene gene = Gene(gene_length,&cities,true);
    Gene new_gene = Gene(gene_length,&cities,true);
    Gene best_gene = gene;

    output.open("evolution.dat");
    first.open("first.gene");
    for (uint i = 0; i < gene.GetLength(); i++){    
        City c = gene.GetElement(i);
        first << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;
    }
    first.close();
    
    for (int i = temperature_steps; i > 0; i--){
        temp = (temp_init/(double)temperature_steps)*i;
        if(i%100==0)cout<< "---------- Temperature " << temp << " ----------" << endl;
        accepted_1 = 0;accepted_2 = 0;accepted_3 = 0;
        for (int j = 0; j < steps_per_T; j++) {
            new_gene = gene;
            new_gene.Mutation_Pair();
            A = min(1.,exp(-(new_gene.Cost()-gene.Cost())/temp));
            if(Random()<A){ gene = new_gene; accepted_1 ++;}
            if(best_gene.Cost()>gene.Cost())best_gene=gene;

            new_gene = gene;
            new_gene.Mutation_Inversion();
            A = min(1.,exp(-(new_gene.Cost()-gene.Cost())/temp));
            if(Random()<A){ gene = new_gene; accepted_2 ++;}
            if(best_gene.Cost()>gene.Cost())best_gene=gene;
            /*
            new_gene = gene;
            new_gene.Mutation_SwitchEnd();
            A = min(1.,exp(-(new_gene.Cost()-gene.Cost())/temp));
            if(Random()<A){ gene = new_gene; accepted_3 ++;}
            if(best_gene.Cost()>gene.Cost())best_gene=gene;
            */
        }
        output << setw(wd) << temperature_steps-i <<  setw(wd) << best_gene.GetCost() << endl;
        if(i%100==0){
            cout<< "Pair Acceptation rate at temp " << temp << " is: " << (double)accepted_1/(double)steps_per_T << endl; 
            cout<< "Inversion Acceptation rate at temp " << temp << " is: " << (double)accepted_2/(double)steps_per_T << endl; 
            //cout<< "End Acceptation rate at temp " << temp << " is: " << (double)accepted_3/(double)steps_per_T << endl << endl; 

            cout<< "Best gene at temp " << temp << " has cost: " << best_gene.Cost() << endl << endl; 
        }
    }
    cout<< "**** Best gene has cost: " << best_gene.Cost() << endl << endl; 
    best.open("best.gene");
    for (uint i = 0; i < best_gene.GetLength(); i++){    
        City c = best_gene.GetElement(i);
        best << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;
    }
    best.close();
    return 0;
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