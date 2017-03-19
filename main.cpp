#include "mpi.h"
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

vactor<double> findVectorOfCoefB(vector<int> &, vector<double> &, vector<double> &);
vector<vector<double>> transpose(vector<vector<double>> &);
vector<vector<double>> multiply(vector<vector<double>> &, vector<vector<double>> &);

vector<int> DEFAULT_COEFF_B(7, 1);
vector<double> calculate(const vector<int>&, double, const vector<int>& = DEFAULT_COEFF_B);

int main(int argc, char * argv[]) {
    int numtasks, taskid;

    MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks); 
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

    vector<double> yA;
    vector<double> xA;
    vector<double> yB;
    vector<double> xB;
    
    ifstream ifile("data.txt");

    if (!ifile.is_open()) {
        cerr << "There was a problem opening the input file!\n";
        exit(1);
    }

    double t = 0.0;
    double indxCounter = 0.0;
    while(ifile >> t) {
        indxCounter++;
        if(yA.size() > yB.size()) {
            yA.push_back(t);
            xA.push_back(indxCounter);
        } else {
            yB.push_back(t);
            xB.push_back(indxCounter);
        }
    }

   vector<int> s(5,4);
    //cout << firstT.size() << secondT.size()<< endl;
    calculate(s, 2);
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

vector<double> calculate(const vector<int>& indxFunc, double x, const vector<int>& b) {

    vector<double> calculatedVal;

    for(int i = 0; i < indxFunc.size(); i++) {
        switch(indxFunc[i]) {
            case 1:
                calculatedVal.push_back(b[0] * 1);
            break;
            
            case 2:
                calculatedVal.push_back(b[1] * sin(4 * x));
            break;

            case 3:
               calculatedVal.push_back(b[2] * tan(7 * x));
            break;

            case 4:
                 calculatedVal.push_back(b[3] * cos(12 * x));
            break;

            case 5:
                 calculatedVal.push_back(b[4] * (1 / tan(7 * x)));
            break;

            case 6:
                calculatedVal.push_back(b[5] * sin(7 * x));
            break;

            case 7:
                calculatedVal.push_back(b[6] * cos(4 * x));
            break;

            default:
            break;
        }
    }

    // for(int i = 0; i < calculatedVal.size(); i++)
    //     cout << calculatedVal[i] << " ";

    // cout << endl;

    return calculatedVal;
}