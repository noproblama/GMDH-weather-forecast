#include "mpi.h"
#include <fstream>
#include <vector>
#include <iostream>

int main(int argc, char * argv[]) {
    int numtasks, taskid;

    MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks); 
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

    std::vector<double> firstT;
    std::vector<double> secondT;

    std::ifstream ifile("data.txt");

    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);
    }

    double t = 0.0;
    while(ifile >> t) {
        if(firstT.size() > secondT.size())
            secondT.push_back(t);
        else 
            firstT.push_back(t);
    }

   
    std::cout << firstT.size() << secondT.size()<< std::endl;

    MPI_Finalize();
}