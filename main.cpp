#include "mpi.h"
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <numeric>
using namespace std;

vector<double> DEFAULT_COEFF_B(7, 1);

vector<double> findVectorOfCoefB(const vector<int> &, const vector<double> &, const vector<double> &);
vector<vector<double> > transpose(const vector<vector<double> > &);
vector<vector<double> > operator*(const vector<vector<double> >& a, const vector<vector<double> >& b);
vector<double> calculate(const vector<int>&, double, const vector<double>& = DEFAULT_COEFF_B);

double CalcDeterminant(vector<vector<double> > , int);
vector<vector<double> > GetMinor(vector<vector<double> >, int, int, int);
vector<vector<double> > inverse(vector<vector<double> >);
string binary(unsigned);
double squareComparation(const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&);


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
    int indxCounter = 0;
    while(ifile >> t) {
        indxCounter++;
        if(yA.size() < yB.size()) {
            yA.push_back(t);
            xA.push_back((double)indxCounter);
        } else {
            yB.push_back(t);
            xB.push_back((double)indxCounter);
        }
    }
            // cout << yA.size() << " ";

            // cout << yB.size() << " ";


    // for (int i = 0; i < yB.size(); ++i)
    // {
    //     cout << xB[i] << " " << yB[i] << endl ;
    // }


   // vector<vector<double> > r1(5, vector<double>(3, 1));
   // vector<vector<double> > r2(3, vector<double>(7, 1));
   // vector<vector<double> > r3(7, vector<double>(3, 1));
   //(r1*r2)*r3;
   
    //cout << firstT.size() << secondT.size()<< endl;

    // double text_xa[] = {1.0, 1.2, 1.3, 1.4, 1.5};
    // double text_ya[] = {1.1, 1.2, 1.3, 1.4, 1.5};

    // vector<double> x(text_xa, text_xa + sizeof(text_xa) / sizeof(double));
    // vector<double> y(text_ya, text_ya + sizeof(text_ya) / sizeof(double));




    vector <vector<int> > L;
    for (int i = 1; i < 64; ++i)
    {
        L.push_back(vector<int>(1,1));
        string binary_number = binary(i);
        for (int j = 0; j < binary_number.size(); ++j) {
            if (binary_number.compare(j,1,"1") == 0) {
                L[i-1].push_back(j + 2);
            }
        }
    }

    vector<vector<double> > Bs;
    for (int i = 0; i < L.size(); ++i)
    {
        Bs.push_back(vector<double> (findVectorOfCoefB(L[i], xA, yA)));  

        // for (int j = 0; j < Bs.size(); ++j)
        // {
        //     cout << Bs[i][j] << " ";

        // }
        // cout << endl; 
    }

    vector<double> square_comparations;
    for (int i = 0; i < Bs.size(); ++i)
    {
        square_comparations.push_back(squareComparation(L[i], Bs[i], xB, yB));
    }
    
    // for (int i = 0; i < square_comparations.size(); ++i)
    // {
    //     cout << square_comparations[i] << " ";
    // }

    // cout << endl;

    int min_pos = distance(square_comparations.begin(),min_element(square_comparations.begin(),square_comparations.end()));
    cout << "The best model is the model #" << min_pos + 1 << " (" << square_comparations[min_pos] << ") "<<" with b: ";
    for (int j = 0; j < Bs[min_pos].size(); ++j)
    {
        cout << Bs[min_pos][j] << " ";

    }
    cout << endl; 


    // for (int j = 0; j < b.size(); ++j)
    // {
    //     cout << b[j] << " ";

    // }
    // cout << endl;


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

    // cout << "COUNTED ONE" << endl;
    // for(int i = 0; i < X_functions.size(); i++){
    //     for (int j = 0; j < X_functions[0].size(); ++j)
    //     {
    //         cout << X_functions[i][j] << " ";
    
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    X_transposed = transpose(X_functions);

    // cout << "TRANSPOSED ONE" << endl;
    // for(int i = 0; i < X_transposed.size(); i++){
    //     for (int j = 0; j < X_transposed[0].size(); ++j)
    //     {
    //         cout << X_transposed[i][j] << " ";
    
    //     }
    //     cout << endl;
    // }

    Y = transpose(YT);

    // vector<vector<double> > hui, pizda, pizda_I, upyachka, upyachka_T;
    // hui = X_transposed*X_functions;
    // pizda = hui*X_transposed;
    // pizda_I = inverse(pizda);
    // upyachka = pizda_I*X_transposed;
    // upyachka_T = transpose(upyachka);

    // cout << "SIZE X: " << X_functions.size() << "x" << X_functions[0],size() << endl;
    // cout << "SIZE XT: " << X_transposed.size() << "x" << X_transposed[0],size() << endl;
    // cout << "SIZE (X*XT)^-1: " << hui.size() << "x" << hui[0],size() << endl;
    // cout << "SIZE X*XT: " << hui.size() << "x" << hui[0],size() << endl;
    vector<double> calculated_bs (transpose((inverse(X_transposed*X_functions)*X_transposed)*Y)[0]);

    // cout<< ((X_transposed*X_functions)*X_transposed) << endl;
    // cout<< Y[0].size() << endl;

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

vector<double> calculate(const vector<int>& indxFunc, double x, const vector<double>& b) {

    vector<double> calculatedVal;

    for(int i = 0; i < indxFunc.size(); i++) {
        switch(indxFunc[i]) {
            case 1:
                calculatedVal.push_back(b[0] * 1);
            break;
            
            case 2:
                calculatedVal.push_back(b[1] * sin(2 * 3.14 * x / 365));
            break;

            case 3:
               calculatedVal.push_back(b[2] * cos(2 * 3.14 * x / 365));
            break;

            case 4:
                 calculatedVal.push_back(b[3] * sin(24 * 3.14 * x / 365));
            break;

            case 5:
                 calculatedVal.push_back(b[4] * cos(24 * 3.14 * x / 365));
            break;

            case 6:
                calculatedVal.push_back(b[5] * sin(12 * 3.14 * x / 365));
            break;

            case 7:
                calculatedVal.push_back(b[6] * cos(12 * 3.14 * x / 365));
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
vector<vector<double> > inverse(vector<vector<double> > matrix)
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

// Convert int to binary string
string binary(unsigned x)
{
    string s;
    do
    {
        s.push_back('0' + (x & 1));
    } while (x >>= 1);

    reverse(s.begin(), s.end());
    for (int i = s.size(); i < 6; i++) {
        s.insert(s.begin(),'0');
       
    }
    return s;

}

double squareComparation(const vector<int> &L, const vector<double> &b, const vector<double> &xB, const vector<double> &yB) {
    double  result,
            model_minus_table_sqr_sum,
            table_sqr_sum;

    for (int i = 0; i < xB.size(); ++i)
    {
        vector<double> model_res_vector = calculate(L, xB[i], b);
        double  model_res = accumulate(model_res_vector.begin(), model_res_vector.end(), 0);
        double  table_res = yB[i];

        model_minus_table_sqr_sum += pow(model_res-table_res, 2);
        table_sqr_sum += pow(table_res, 2);
    }

    result = model_minus_table_sqr_sum / table_sqr_sum;
    
    return result;
}