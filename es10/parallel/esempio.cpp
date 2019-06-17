#include "mpi.h"
#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
    MPI::Init(argc,argv);
    int size = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();
    if(size>3){
        cout<<"Hai scelto troppi processi"<<endl;
        return 1;}
    int irecv[3];
    for(int i=0;i<3;i++) irecv[i]=0;
    int isend = rank + 1;
    MPI_Gather(&isend,1,MPI_INTEGER,irecv,1,MPI_INTEGER,0,MPI::COMM_WORLD);
    if(rank==0) cout<< "irecv: " <<irecv[0] <<" "<<irecv[1] <<" " <<irecv[2] <<endl;
    MPI::Finalize();
    return 0;
}