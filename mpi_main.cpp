#include "mpi.h"
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <numeric>
using namespace std;

const vector<double> DEFAULT_COEFF_B(7, 1);
const double PI = 3.14159265359;

vector<double> findVectorOfCoefB(const vector<int> &, const vector<double> &, const vector<double> &);
vector<vector<double> > transpose(const vector<vector<double> > &);
vector<vector<double> > operator*(const vector<vector<double> >& a, const vector<vector<double> >& b);
vector<double> calculate(const vector<int>&, double, const vector<double>& = DEFAULT_COEFF_B);
double CalcDeterminant(vector<vector<double> >, int);
vector<vector<double> > GetMinor(vector<vector<double> >, int, int, int);
vector<vector<double> > inverse(vector<vector<double> >);
string binary(unsigned);
double squareComparation(const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&);
double differenceComparation(const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&);

int main(int argc, char * argv[]) {
    int numtasks, taskid;

    MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks); 
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    


    //cout << xTable[12] << " " << yTable[12] << endl;
    /*<<<read from file*/

    /*>>>generate all variants of models*/
    // for (int i = 1; i < 64; ++i)
    // {
    //     L.push_back(vector<int>(1,1));
    //     string binary_number = binary(i);
    //     for (int j = 0; j < binary_number.size(); ++j) {
    //         if (binary_number.compare(j,1,"1") == 0) {
    //             L[i-1].push_back(j + 2);
    //         }
    //     }
    // }
    /*<<<generate all variants of models*/

    
    if(taskid != 0) {

        vector<double> yTable;
        vector<double> xTable;
        vector<double> yA;
        vector<double> xA;
        vector<double> yB;
        vector<double> xB;
        vector <vector<int> > L;

        /*>>>read from file*/
        ifstream ifile("test2.txt");

        if (!ifile.is_open()) {
            cerr << "There was a problem opening the input file!\n";
            exit(0);
        }
cout << "read data" << endl;
        double t = 0.0;
        int indxCounter = 0;
        while(ifile >> t) {
            indxCounter++;
            if(yA.size() < yB.size()) {
                yA.push_back(t);
                xA.push_back(indxCounter);
            } else {
                yB.push_back(t);
                xB.push_back(indxCounter);
            }
            xTable.push_back(indxCounter);
            yTable.push_back(t);
        }

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

        int nModel = L.size() / (numtasks - 1);
        cout << "nModel: " << nModel << endl;

        int start = (taskid - 1) * nModel;
        int end = taskid * nModel;
        cout << "start: " << start << " end: " << end << endl;
        vector<vector<double> > bA;
        
        for (int i = start; i < end; i++)
            bA.push_back(vector<double> (findVectorOfCoefB(L[i], xA, yA)));

        vector<vector<double> > bB;

        for (int i = start; i < end; i++)
            bB.push_back(vector<double> (findVectorOfCoefB(L[i], xB, yB)));

        vector<double> square_comparations;
        for (int i = start; i < end; i++){
            square_comparations.push_back(squareComparation(L[i], bA[i - start], xB, yB)); //maybe fail
            //cout << "i - start: " << i - start << endl;
        }
       
        int min_pos_sqr = distance(square_comparations.begin(), min_element(square_comparations.begin(), square_comparations.end()));

        MPI_Send(&min_pos_sqr, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&bA[min_pos_sqr][0], L[min_pos_sqr].size(), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);


        vector<double> bias_comparations;
        if(bA.size() != bB.size()) {
            cout << "bA and bB must have the same length";
            exit(0);
        }

        for (int i = start; i < end; i++)
            bias_comparations.push_back(differenceComparation(L[i], bA[i - start], bB[i - start], xTable, yTable));

        int min_pos_bias = distance(bias_comparations.begin(), min_element(bias_comparations.begin(), bias_comparations.end()));
        
        MPI_Send(&min_pos_bias, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
        MPI_Send(&bA[min_pos_bias][0], L[min_pos_bias].size(), MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
        MPI_Send(&bB[min_pos_bias][0], L[min_pos_bias].size(), MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
    }

    if(taskid == 0) {
        vector<double> yTable;
        vector<double> xTable;
        vector<double> yA;
        vector<double> xA;
        vector<double> yB;
        vector<double> xB;
        vector <vector<int> > L;

        /*>>>read from file*/
        ifstream ifile("test3.txt");

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
                xA.push_back(indxCounter);
            } else {
                yB.push_back(t);
                xB.push_back(indxCounter);
            }
            xTable.push_back(indxCounter);
            yTable.push_back(t);
        }

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



        vector<vector<double> > results_sqr;
        vector<vector<double> > results_biasA;
        vector<vector<double> > results_biasB;
        vector<int> index_sqr;
        vector<int> index_bias;

        for(int p = 1; p < numtasks; p++) {
            int min_pos_sqr, min_pos_bias;
            
            MPI_Recv(&min_pos_sqr, 1, MPI_INT, p, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            vector<double> square(L[min_pos_sqr].size());
            index_sqr.push_back(min_pos_sqr);
            MPI_Recv(&square[0], L[min_pos_sqr].size(), MPI_DOUBLE, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            results_sqr.push_back(square);

            MPI_Recv(&min_pos_bias, 1, MPI_INT, p, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            index_bias.push_back(min_pos_bias);
            vector<double> biasA(L[min_pos_bias].size());
            vector<double> biasB(L[min_pos_bias].size());
            MPI_Recv(&biasA[0], L[min_pos_bias].size(), MPI_DOUBLE, p, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            results_biasA.push_back(biasA);
            MPI_Recv(&biasB[0], L[min_pos_bias].size(), MPI_DOUBLE, p, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            results_biasB.push_back(biasB);                    
        }

        int nModelTail = L.size() - ((L.size() / (numtasks - 1)) * (numtasks - 1));
        cout << "MASTER nModel: " << nModelTail << endl;
        if(nModelTail > 0) {
            // int start = L.size() - nModelTail - 1;
            // int end = L.size();

            // vector<vector<double> > bA;
            // for (int i = start; i < end; i++)
            //     bA.push_back(vector<double> (findVectorOfCoefB(L[i], xA, yA)));

            // vector<vector<double> > bB;
            // for (int i = start; i < end; i++)
            //     bB.push_back(vector<double> (findVectorOfCoefB(L[i], xB, yB)));

            // vector<double> square_comparations;
            // for (int i = start; i < end; i++)
            //     square_comparations.push_back(squareComparation(L[i], bA[i - start], xB, yB)); //maybe fail
        
            // int min_pos_sqr = distance(square_comparations.begin(), min_element(square_comparations.begin(), square_comparations.end()));

            // results_sqr.push_back(bA[min_pos_sqr]);
            // index_sqr.push_back(min_pos_sqr);

            // vector<double> difference_comparations;
            // if(bA.size() != bB.size()) {
            //     cout << "bA and bB must be the same length";
            //     exit(0);
            // }

            // for (int i = start; i < end; i++)
            //     difference_comparations.push_back(differenceComparation(L[i], bA[i - start], bB[i - start], xTable, yTable));

            // int min_pos_bias = distance(difference_comparations.begin(), min_element(difference_comparations.begin(), difference_comparations.end()));
            // results_biasA.push_back(bA[min_pos_bias]);
            // results_biasB.push_back(bB[min_pos_bias]);
            // index_bias.push_back(min_pos_bias);
        }

        /*Final comparation*/
        vector<double> square_final;
        for (int i = 0; i < results_sqr.size(); i++)
            square_final.push_back(squareComparation(L[index_sqr[i]], results_sqr[i], xB, yB)); //maybe fail
       
        int min_pos_sqr_f = distance(square_final.begin(), min_element(square_final.begin(), square_final.end()));

        cout << "The best model by SQR criteria is the model #" << index_sqr[min_pos_sqr_f] + 1 << " (" << square_final[min_pos_sqr_f] << ") "<<" with b: ";
        for (int j = 0; j < results_sqr[min_pos_sqr_f].size(); ++j)
            cout << results_sqr[min_pos_sqr_f][j] << " ";
        cout << "// with model functions #: ";

        for (int j = 0; j < L[index_sqr[min_pos_sqr_f]].size(); ++j)
            cout << L[index_sqr[min_pos_sqr_f]][j] << " ";
        cout << endl;  
        
        vector<double> bias_final;
      
        if(results_biasA.size() != results_biasB.size()) {
            cout << "Troubles" << endl;
            exit(0);
        }
        for (int i = 0; i < results_biasA.size(); i++)
            bias_final.push_back(differenceComparation(L[index_bias[i]], results_biasA[i], results_biasB[i], xTable, yTable));

        int min_pos_bias_f = distance(bias_final.begin(), min_element(bias_final.begin(), bias_final.end()));

        cout << "The best model by BIAS criteria is the model #" << index_bias[min_pos_bias_f] + 1 << " (" << bias_final[min_pos_bias_f] << ") "<<" with b: ";
        for (int j = 0; j < results_biasA[min_pos_bias_f].size(); ++j)
            cout << results_biasA[min_pos_bias_f][j] << " ";
        cout << "// with model functions #: ";

        for (int j = 0; j < L[index_bias[min_pos_bias_f]].size(); ++j)
            cout << L[index_bias[min_pos_bias_f]][j] << " ";
        cout << endl;  
        

    }

    
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
    vector<double> calculated_bs (transpose((inverse((X_transposed*X_functions))*X_transposed)*Y)[0]);

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
    vector<double> calcVal;

    for(int i = 0; i < indxFunc.size(); i++) {
        switch(indxFunc[i]) {
            case 1:
                calcVal.push_back(b[0] * 1);
            break;
            
            case 2:
                calcVal.push_back(b[1] * sin(2 * PI * x / 365));
            break;

            case 3:
               calcVal.push_back(b[2] * cos(2 * PI * x / 365));
            break;

            case 4:
                 calcVal.push_back(b[3] * 5 * x);
            break;

            case 5:
                 calcVal.push_back(b[4] * cos(24 * PI * x / 365));
            break;

            case 6:
                calcVal.push_back(b[5] * sin(PI * x / 14));
            break;

            case 7:
                calcVal.push_back(b[6] * cos(PI * x / 14));
            break;

            default:
            break;
        }
    }

    // for(int i = 0; i < calcVal.size(); i++)
    //     cout << calcVal[i] << " ";

    // cout << endl;

    return calcVal;
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
            model_minus_table_sqr_sum = 0.0,
            table_sqr_sum =0.0;

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

double differenceComparation(const vector<int> &L, const vector<double> &bA, const vector<double> &bB, const vector<double> &xTable, const vector<double> &yTable) {
    double  result,
            models_diff_sqr_sum = 0.0,
            table_sqr_sum = 0.0;

    for (int i = 0; i < xTable.size(); ++i)
    {
        vector<double> model_res_vector_A = calculate(L, xTable[i], bA);
        vector<double> model_res_vector_B = calculate(L, xTable[i], bB);
        double  model_res_A = accumulate(model_res_vector_A.begin(), model_res_vector_A.end(), 0);
        double  model_res_B = accumulate(model_res_vector_B.begin(), model_res_vector_B.end(), 0);

        models_diff_sqr_sum += pow(model_res_A-model_res_B, 2);

        double  table_res = yTable[i];    
        table_sqr_sum += pow(table_res, 2);
    }

    result = models_diff_sqr_sum / table_sqr_sum;
    
    return result;
}