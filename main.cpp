#include "mpi.h"
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
using namespace std;

vector<int> DEFAULT_COEFF_B(7, 1);

vector<double> findVectorOfCoefB(const vector<int> &, const vector<double> &, const vector<double> &);
vector<vector<double> > transpose(const vector<vector<double> > &);
vector<vector<double> > operator*(const vector<vector<double> >& a, const vector<vector<double> >& b);
vector<double> calculate(const vector<int>&, double, const vector<int>& = DEFAULT_COEFF_B);

double CalcDeterminant(vector<vector<double> > , int);
vector<vector<double> > GetMinor(vector<vector<double> >, int, int, int);
vector<vector<double> > MatrixInversion(vector<vector<double> >);

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
        exit(0);
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

   vector<vector<double> > r1(5, vector<double>(3, 1));
   vector<vector<double> > r2(3, vector<double>(7, 1));
   vector<vector<double> > r3(7, vector<double>(3, 1));
   //(r1*r2)*r3;
   
    //cout << firstT.size() << secondT.size()<< endl;

    int test_l[] = {1, 6};
    double text_xa[] = {1.0, 1.2, 1.3, 1.4, 1.5};
    double text_ya[] = {1.1, 1.2, 1.3, 1.4, 1.5};

    vector<int> L(test_l, test_l + sizeof(test_l) / sizeof(int));
    vector<double> x(text_xa, text_xa + sizeof(text_xa) / sizeof(double));
    vector<double> y(text_ya, text_ya + sizeof(text_ya) / sizeof(double));


    vector<double> b (findVectorOfCoefB(L, x, y));

    MPI_Finalize();
}

// Find coefficients of b for the model
// L - indexes of used functions
// XA - 1d vector of x from A
// YA - 1d vector of y from A
vector<double> findVectorOfCoefB(const vector<int> &L, const vector<double> &XA, const vector<double> &YA) {
    vector<vector<double> >     X_functions,    // matrix of the model's functions results
                                X_transposed,
                                YT(1, YA),
                                Y;
    for (int i = 0; i < XA.size(); ++i)
    {
        X_functions.push_back(calculate(L, XA[i]));
    }                    

    cout << "COUNTED ONE" << endl;
    for(int i = 0; i < X_functions.size(); i++){
        for (int j = 0; j < X_functions[0].size(); ++j)
        {
            cout << X_functions[i][j] << " ";
    
        }
        cout << endl;
    }
    cout << endl;

    X_transposed = transpose(X_functions);

    cout << "TRANSPOSED ONE" << endl;
    for(int i = 0; i < X_transposed.size(); i++){
        for (int j = 0; j < X_transposed[0].size(); ++j)
        {
            cout << X_transposed[i][j] << " ";
    
        }
        cout << endl;
    }

    Y = transpose(YT);
    vector<double> calculated_bs (transpose((((X_transposed*X_functions)*X_transposed)*Y))[0]);

    // cout<< ((X_transposed*X_functions)*X_transposed) << endl;
    // cout<< Y[0].size() << endl;

    // cout << "bs calculated:" << endl;

    // for (int j = 0; j < calculated_bs.size(); ++j)
    // {
    //     cout << calculated_bs[j] << " ";

    // }
    // cout << endl;

    return calculated_bs;
}

// Transpose a 2d matrix
vector<vector<double> > transpose(const vector<vector<double> > &matrix) {
    int size[] = {(int)matrix.size(), (int)matrix[0].size()};
    vector<vector<double> > transposedMatrix(size[1], vector<double> (size[0], 0));
    for (int i = 0; i < matrix.size(); ++i)
    {
        for (int j = 0; j < matrix[i].size(); ++j)
        {
            transposedMatrix[j][i] = matrix[i][j];
        }
    }
    return transposedMatrix;
};

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

vector<vector<double> > operator*(const vector<vector<double> >& a, const vector<vector<double> >& b)
{
    vector<vector<double> > result(a.size(), vector<double>(b[0].size()));
    for(size_t i = 0; i < a.size(); i++)
    {
        for(size_t j = 0; j < b[0].size(); j++)
        {
            double sum = 0;
            for(size_t k = 0; k < a[0].size(); k++)
            {
                sum += a[i][k] * b[k][j];
            }
            result[i][j] = sum;
        }
    }

    // for(size_t i = 0; i < result.size(); i++)
    // {
    //     for(size_t j = 0; j < result[0].size(); j++)
    //     {
    //         cout << result[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // cout << endl;

    return  result;
}

// matrix inversioon
vector<vector<double> > MatrixInversion(vector<vector<double> > matrix)
{
    
    int order = matrix.size();
    vector<vector<double> > res(order, vector<double>(order));
    // get the determinant of matrix
    double det = 1.0 / CalcDeterminant(matrix, order);
 
    // memory allocation
    vector<vector<double> > minorM(order-1, vector<double>(order-1));

 
    for(int j=0;j<order;j++)
    {
        for(int i=0;i<order;i++)
        {
            // get the co-factor (matrix) of A(j,i)
            minorM = GetMinor(matrix,j,i,order);
            res[i][j] = det*CalcDeterminant(minorM,order-1);
            if( (i+j)%2 == 1)
                res[i][j] = -res[i][j];
        }
    }

    return res;
}
 
// calculate the cofactor of element (row,col)
vector<vector<double> > GetMinor(vector<vector<double> > src, int row, int col, int order)
{
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;
    vector<vector<double> > res(order, vector<double>(order));
 
    for(int i = 0; i < order; i++ )
    {
        if( i != row )
        {
            colCount = 0;
            for(int j = 0; j < order; j++ )
            {
                // when j is not the element
                if( j != col )
                {
                    res[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }
 
    return res;
}
 
// Calculate the determinant recursively.
double CalcDeterminant(vector<vector<double> > matrix, int order)
{
    //order must be >= 0
    //stop the recursion when matrix is a single element
    if( order == 1 )
        return matrix[0][0];
 
    // the determinant value
    double det = 0;
 
    // allocate the cofactor matrix
    vector<vector<double> > minorM(order-1, vector<double>(order - 1));
 
    for(int i = 0; i < order; i++ )
    {
        // get minor of element (0,i)
        minorM = GetMinor( matrix, 0, i , order);
        // the recusion is here!
 
        det += (i%2==1?-1.0:1.0) * matrix[0][i] * CalcDeterminant(minorM,order-1);
        
    }

 
    return det;
}