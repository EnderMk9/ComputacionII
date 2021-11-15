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
void coutmat(Matrix& M){                  // Input is Matrix M passed with & as a pointer
    cout << endl;
    int rows=M.size(); int cols=M[0].size(); // size of M
    for (int i = 0; i < rows; i++){          // for every row
        for (int j = 0; j < cols; j++){      // for every column
            cout << M[i][j] << " ";          // Show Mij index of M followed by a space
        }cout << endl;                       // next line
    } cout << endl;
}

// displays a vector in console
void coutvec(Vector& v){                  // Input is Matrix M passed with & as a pointer
    int rows=v.size();                    // Size of the vector.
    for (int i = 0; i < rows; i++){       // for every row
            cout << v[i] << " ";          // Show Mij index of M followed by a space
    } cout << endl;
}

void Scoutvec(SVector& v){                  // Input is Matrix M passed with & as a pointer
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
