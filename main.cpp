#include "mpi.h"
#include <fstream>
#include <vector>
#include <iostream>
using namespace std;

vactor<double> findVectorOfCoefB(vector<int> &, vector<double> &, vector<double> &);
vector<vector<double>> transpose(vector<vector<double>> &);
vector<vector<double>> multiply(vector<vector<double>> &, vector<vector<double>> &);


int main(int argc, char * argv[]) {
    int numtasks, taskid;

    MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks); 
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

    vector<double> firstT;
    vector<double> secondT;

    ifstream ifile("data.txt");

    if (!ifile.is_open()) {
        cerr << "There was a problem opening the input file!\n";
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

// Find coefficients of b for the model
// L - indexes of used functions
// XA - 1d vector of x from A
// YA - 1d vector of y from A
vactor<double> findVectorOfCoefB(vector<int> &L, vector<double> &XA, vector<double> &YA) {
    vector<double>  calculated_bs, 
                    X_functions,    // matrix of the model's functions results
                    X_transposed; 
    for (vector<double>::iterator x = XA.begin(); i != XA.end(); ++x)
    {   
        X_functions.push_back(calculate(L, *x));
    }

    X_transposed = transpose(X_functions);

    calculated_bs = multiply(multiply(multiply(X_transposed, X_functions),X_transposed), YA);

    return calculated_bs;
}

// Transpose 2d matrix
vector<vector<double>> transpose(vector<vector<double>> &matrix) {
    int size[] = {matrix.size(), matrix[0].size()};
    vector<vector<double>> transponedMatrix(size[1], vector<double> (size[0], 0));
    for (int i = 0; i < matrix.size(); ++i)
    {
        for (int j = 0; j < matrix[i].size(); ++j)
        {
            transponedMatrix[j][i] = matrix[i][j];
        }
    }
    return transponedMatrix;
};

vector<vector<double>> multiply(vector<vector<double>> &matrixA, vector<vector<double>> &matrixB) {
    vector<vector<double>> multipliedM;
    return multipliedM;
}