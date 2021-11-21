#include <vector>
#include <string>

typedef vector<vector<double>> Matrix;    // Create a type called Matrix that is a vector of vectors
typedef vector<double> Vector;            // Create a type called Vector that is a vector

// The next two are the same but the entries are int instead of double
typedef vector<vector<int>> IMatrix;
typedef vector<int> IVector;

// The next two are the same but the entries are long double instead of double
typedef vector<vector<long double>> LMatrix;
typedef vector<long double> LVector;

// The next two are the same but the entries are long double instead of double
typedef vector<vector<string>> SMatrix;
typedef vector<string> SVector;

//------------------------------------------
// Cout vectors
//------------------------------------------

// displays a matrix in console
void coutmat(Matrix& M, int pres = 16){      // Input is Matrix M passed with & as a pointer
    cout.precision(pres);
    cout << endl;
    int rows=M.size(); int cols=M[0].size(); // size of M
    for (int i = 0; i < rows; i++){          // for every row
        for (int j = 0; j < cols; j++){      // for every column
            cout << M[i][j] << " ";          // Show Mij index of M followed by a space
        }cout << endl;                       // next line
    } cout << endl;
}

// displays a vector in console
void coutvec(Vector& v, int pres = 16){   // Input is Matrix M passed with & as a pointer
    int rows=v.size();                    // Size of the vector.
    cout.precision(pres);
    for (int i = 0; i < rows; i++){       // for every row
            cout << v[i] << " ";          // Show Mij index of M followed by a space
    } cout << endl;
}

// displays a vector in console
void Icoutvec(IVector& v){   // Input is Matrix M passed with & as a pointer
    int rows=v.size();                    // Size of the vector.
    for (int i = 0; i < rows; i++){       // for every row
            cout << v[i] << " ";          // Show Mij index of M followed by a space
    } cout << endl;
}

void Scoutvec(SVector& v){                // Input is Matrix M passed with & as a pointer
    int rows=v.size();                    // Size of the vector.
    for (int i = 0; i < rows; i++){       // for every row
            cout << v[i] << " ";          // Show Mij index of M followed by a space
    } cout << endl;
}

//------------------------------------------
// Creation and conversion tools
//------------------------------------------

// Creates an n-vector with all components with a certain value
Vector VecFull(double value, int n){
    Vector v(n,value);
    return v;
}

// Creates a matrix of n rows and m columns with all components having the same value
Matrix MatFull(double value, int n, int m){
    Matrix M(n,Vector(m,value));
    return M;
}

//------------------------------------------
// Elementary operations
//------------------------------------------

// transposes a square matrix
Matrix transposeSqr(Matrix& M){ // Input is a square Matrix M passed with & as a pointer
    int rows=M.size(); int cols=M[0].size();  // Size of M
    if (cols != rows){                        // Check if M is square
        cout << "NOT SQUARE" << endl;         // error
        return Matrix {0};                    // break
        }
    for (int i = 0; i < rows; i++){           // for every row
        for (int j = 0; j <= i; j++){         // for every column
            swap(M[i][j],M[j][i]);            // swap the entries
        }
    }
    return M;
}

// transposes a matrix
Matrix transpose(Matrix& M){ // Input is a square Matrix M passed with & as a pointer
    int rows=M.size(); int cols=M[0].size();  // Size of M
    Matrix T( cols,vector<double>(rows,0));   // Initialize the transpose
    for (int i = 0; i < cols; i++){           // for every row
        for (int j = 0; j < rows; j++){       // for every column
            T[i][j] = M[j][i];                // Define the transpose
        }
    }
    return T;
}

// Multiplies a vector A by a scalar lambda
Vector ScalMult(Vector& A, double lambda){
    int n = A.size();
    Vector B(n,0);                          // Define the solution
    for (int i = 0; i < n; i++){            // For every component
        B[i] = lambda*A[i];
    }
    return B;
}

// Returns the component-wise sum of two vectors
// A and B must be the same size
Vector VecSum(Vector& A, Vector& B){   // A_i*B_i
    int An = A.size(); int Bn = B.size();    // Size of each vector
    if (An != Bn){                           // Check if they are the same size
        cout << "NOT SAME SIZE" << endl;     // error
        return {};}                          // break
    Vector C(An,0);                          // Define the solution
    for (int i = 0; i < An; i++){            // For every component
        C[i] = A[i]+B[i];                    // Multiplication
    }
    return C;
}

// Returns the difference of two vectors
// A and B must be the same size
Vector VecDiff(Vector& A, Vector& B){  // A-B
    int An = A.size(); int Bn = B.size();    // Size of each vector
    if (An != Bn){                           // Check if they are the same size
        cout << "NOT SAME SIZE" << endl;     // error
        return {};}                          // break
    Vector C(An,0);                          // Define the solution
    for (int i = 0; i < An; i++){            // For every component
        C[i] = A[i]-B[i];                    // Difference
    }
    return C;
}

// Multiplies two vector component-wise
// A and B must be the same size
Vector VecMult(Vector& A, Vector& B){   // A_i*B_i
    int An = A.size(); int Bn = B.size();    // Size of each vector
    if (An != Bn){                           // Check if they are the same size
        cout << "NOT SAME SIZE" << endl;     // error
        return {};}                          // break
    Vector C(An,0);                          // Define the solution
    for (int i = 0; i < An; i++){            // For every component
        C[i] = A[i]*B[i];                    // Multiplication
    }
    return C;
}

// Returns the component-wise division of two vectors
// A and B must be the same size
Vector VecDiv(Vector& A, Vector& B){  // A_i/B_i
    int An = A.size(); int Bn = B.size();    // Size of each vector
    if (An != Bn){                           // Check if they are the same size
        cout << "NOT SAME SIZE" << endl;     // error
        return {};}                          // break
    Vector C(An,0);                          // Define the solution
    for (int i = 0; i < An; i++){            // For every component
        C[i] = A[i]/B[i];                    // Divission
    }
    return C;
}


double VecProd(Vector& A){
    int n = A.size(); double p = 1;
    for (int i = 0; i < n; i++){
        p = p*A[i];
    }
    return p;
}

double DiagonalProd(Matrix& U){             
    int n=U.size(); double D = U[0][0];
    for (int i = 1; i < n; i++){
        D = D*U[i][i];
    }
    return D;
}