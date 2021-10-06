// requieres <vector> and <math.h>
typedef vector< vector<double> > Matrix;  // Create a type called Matrix that is a vector of vectors
typedef vector<double> Vector;            // Create a type called Vector that is a vector

// displays in console a matrix
void coutmat(Matrix& M){                  // Input is Matrix M passed with & as a pointer
    cout << endl;     
    int rows=M.size(); int cols=M[0].size(); // size of M                  
    for (int i = 0; i < rows; i++){          // for every row
        for (int j = 0; j < cols; j++){      // for every column
            cout << M[i][j] << " ";          // Show Mij index of M followed by a space
        }cout << endl;                       // next line
    } cout << endl;
}
// displays in console a matrix
void coutvec(Vector& v){                  // Input is Matrix M passed with & as a pointer
    cout << endl;                          
    int rows=v.size();                        // Size of the vector.
    for (int i = 0; i < rows; i++){          // for every row
            cout << v[i] << " ";          // Show Mij index of M followed by a space
    } cout << endl;
}

// Multiplies two matrices
Matrix matprod(Matrix& A, Matrix& B){  // Inputs are Matrices A and B passed with & as a pointer
    int Arows=A.size(); int Acols=A[0].size();  // sizes of A
    int Brows=B.size(); int Bcols=B[0].size();  // sizes of B
    if (Acols != Brows){                        // compatible multiplication check
        cout << "WRONG MULTIPLICATION" << endl; // error
        return Matrix {0};                      // break
        }
    int Crows=Arows; int Ccols=Bcols;           // dimensions of the product
    Matrix C( Crows,vector<double>(Ccols,0));   // definition of the product
    for (int i = 0; i < Crows; i++){            // every row
        for (int j = 0; j < Ccols; j++){        // every column                      
            for (int k = 0; k < Crows; k++){    // over the common index
                C[i][j] += A[i][k]*B[k][j];     // multiplication and sum
            }
        }
    }
    return C;
}

// transposes a square matrix
Matrix transpose(Matrix M){ // Input is a square Matrix M passed with & as a pointer
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

// Calculates the standard dot product between two vectors
double dotprod(Vector& X,Vector& Y){        // Inputs are two vectos X and Y of the same size, passed with & as a pointer
    int Xsize=X.size(); int Ysize=Y.size(); // Size of each vector
    if (Xsize != Ysize){                    // Of sizes are not equal
        cout << "NOT SAME SIZE" << endl;    // Error
        return 0;                           // Break
        }
    double prod = 0;                        // Define the output
    for (int i = 0; i < Xsize; i++){        // For every component
        prod += X[i]*Y[i];                  // We multiply the entries and sum
    }
    return prod;
}

// calculates the standard norm of a vector
double norm(Vector& X){         // Input is a Vector, passed with & as a pointer
    return sqrt(dotprod(X,X));  // Calculate the norm as sqrt(X*X)
}

// calculates the angle between two vectors
double angle(Vector& X,Vector& Y){ // Inputs are two vectors X and Y, passed with & as a pointer
    double cos_a = dotprod(X,Y)/(norm(X)*norm(Y)); // definition of the cosine of the angle
    double angle = acos(cos_a);                    // acos to obtain the angle (in radians)
    return angle;
}

// Transforms a matrix on to an orthogonal one (if square),
// orthogonalizing the columns interpreted as vectors
Matrix gramschmidt(Matrix& M){
    // M contains vectors as columns
    int d = M.size();    // Number of dimensions in each vector
    int n = M[0].size(); // Number of vectors
    Matrix G = transpose(M);  // Definition of the normal matrix, transpose to use the rows as vectors.
    for (int i = 1; i < n; i++){     // For each vector beggining by the 2nd
        for (int j = 0; j < i; j++){ // For every previous vector
            double jip; double jsqr;
            jip  = dotprod(G[j],G[i]);
            jsqr = dotprod(G[j],G[j]);
            if (jsqr == 0){
                cout << "ERROR NOT INDEPENDENT" << endl;
                return Matrix {0};
            }      
            double proyection =jip/jsqr; // orthogonal proyection of veci onto vec_j
            for (int k = 0; k < d; k++){        // for every component
                G[i][k] -= G[j][k]*proyection;  // Gram-Schmidt
            }
        }
    }
    for (int i = 0; i < d; i++){
        double norm = sqrt(dotprod(G[i],G[i]));
        if (norm == 0){
            cout << "ERROR NOT INDEPENDENT" << endl;
            return Matrix {0};
        }
        
        for (int j = 0; j < n; j++){
            G[i][j] = G[i][j]/norm;
        }
    }
    G = transpose(G); // to have vectors as columns
    return G;
}

// Matrix L{} and Matrix U{} must be zeroes
// A must be square and non-singular
void LU(Matrix& A, Matrix& L, Matrix& U){
    int n=A.size(); int cols=A[0].size();  // Size of M
    if (cols != n){                        // Check if M is square
        cout << "NOT SQUARE" << endl;      // error
        return;}                           // break
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i==j) L[i][j]=1;
            if (i > j){
                L[i][j]=A[i][j];
                for (int k = 0; k < j; k++){
                    L[i][j]-=L[i][k]*U[k][j];}
                L[i][j] = L[i][j]/U[j][j];
            }else if (i <= j){
                U[i][j]=A[i][j];
                for (int k = 0; k < i; k++){
                    U[i][j]-=L[i][k]*U[k][j];}
            }
        }
    }   
}
Vector LUSolve(Matrix& A, Vector& b){
    int n=A.size(); int cols=A[0].size();  // Size of M
    if (cols != n){                        // Check if M is square
        cout << "NOT SQUARE" << endl;      // error
        return {};}                        // break
    Matrix L( n,vector<double>(n,0));
    Matrix U( n,vector<double>(n,0));
    LU(A,L,U);
    coutmat(L); coutmat(U);
    Vector z(n,0); Vector x(n,0); 
    for(int i = 0; i < n; i++){
        z[i] = b[i];
        for(int j = 0; j < i; j++){
            z[i] -= L[i][j]*z[j];
        }}
    for(int i = 0; i < n; i++){
        x[i] = z[i];
        for(int j = i; j < n; j++){
            x[i] -= U[i][j]*x[j];
            }
        x[i] = x[i]/U[i][i];
        }
    return x;
}