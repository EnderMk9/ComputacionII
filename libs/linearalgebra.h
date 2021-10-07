// requires <vector> and <math.h>
// requires "rwlib.h" ( <fstream> and <string> )
typedef vector<vector<double>> Matrix;  // Create a type called Matrix that is a vector of vectors
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
// displays in console a vector
void coutvec(Vector& v){                  // Input is Matrix M passed with & as a pointer
    int rows=v.size();                    // Size of the vector.
    for (int i = 0; i < rows; i++){       // for every row
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

// Returns the component-wise multiplication of two vectors
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

// calculates the angle between two vectors
double angle(Vector& X,Vector& Y){ // Inputs are two vectors X and Y, passed with & as a pointer
    double cos_a = dotprod(X,Y)/(norm(X)*norm(Y)); // definition of the cosine of the angle
    double angle = acos(cos_a);                    // acos to obtain the angle (in radians)
    return angle;
}

// Transforms a matrix on to an orthogonal one (if square),
// orthogonalizing the columns interpreted as vectors
Matrix gramschmidt(Matrix& M){ // M contains vectors as columns
    int d = M.size();          // Number of dimensions in each vector
    int n = M[0].size();       // Number of vectors
    Matrix G = transpose(M);   // Definition of the normal matrix, transpose to use the rows as vectors.
    for (int i = 1; i < n; i++){     // For each vector beggining by the 2nd
        for (int j = 0; j < i; j++){ // For every previous vector
            double jip; double jsqr;
            jip  = dotprod(G[j],G[i]);
            jsqr = dotprod(G[j],G[j]);
            if (jsqr == 0){
                cout << "ERROR NOT INDEPENDENT" << endl;
                return Matrix {0};
            }
            double proyection =jip/jsqr;        // orthogonal proyection of veci onto vec_j
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

// Decomposes a matrix into upper and lower triangular form
// Matrix L{} and Matrix U{} must be zeroes
// A must be square and non-singular
void LU(Matrix& A, Matrix& L, Matrix& U){
    int n=A.size(); int cols=A[0].size();  // Size of M
    if (cols != n){                        // Check if M is square
        cout << "NOT SQUARE" << endl;      // error
        return;}                           // break
    for (int i = 0; i < n; i++){           // for every row
        for (int j = 0; j < n; j++){       // for every column
            if (i==j) L[i][j]=1;           // if the indices are the same Lower es 1s
            if (i > j){                    // For the indices in L
                L[i][j]=A[i][j];           
                for (int k = 0; k < j; k++){
                    L[i][j]-=L[i][k]*U[k][j];}
                L[i][j] = L[i][j]/U[j][j];
            }else if (i <= j){             // For the indices in U
                U[i][j]=A[i][j];
                for (int k = 0; k < i; k++){
                    U[i][j]-=L[i][k]*U[k][j];}
            }
        }
    }
}
// Solves Ax=b using LU decomposition
// Matrix A must be non-singular, square and b same size as A
Vector LUSolve(Matrix& A, Vector& b){
    int n=A.size(); int cols=A[0].size();  // Size of M
    int bsize = b.size();                  // Size of b
    if (cols != n || bsize != n){                           // Check if M is square
        cout << "NOT SQUARE OR NOT SAME SIZE" << endl;      // error
        return {};}                                         // break
    Matrix L( n,vector<double>(n,0));      // Lower matrix of 0s
    Matrix U( n,vector<double>(n,0));      // Upper matrix of 0s
    LU(A,L,U);                             // LU decomposition
    // coutmat(L); coutmat(U);
    Vector z(n,0); Vector x(n,0);       // Define the vectors to solve
    for(int i = 0; i < n; i++){         // First solve Lz=b
        z[i] = b[i];
        for(int j = 0; j < i; j++){
            z[i] -= L[i][j]*z[j];
        }
      }
    for(int i = n-1; i >= 0; i--){      // Solve Ux=z
        x[i] = z[i];
        for(int j = i+1; j < n; j++){
            x[i] -= U[i][j]*x[j];
            }
        x[i] = x[i]/U[i][i];}
    return x;   // Return the solution
}

// Returns the diagonal of a matrix as a vector
// A must be
Vector MatrixDiag(Matrix& A){
    int n=A.size(); int cols=A[0].size();  // Size of M
    if (cols != n ){                       // Check if M is square
        cout << "NOT SQUARE" << endl;      // error
        return {};}                        // break
    Vector D(n,0);                         // Define the output
    for (int i = 0; i < n; i++){           // for every line or row
        D[i] = A[i][i];}                   // D_i = A_ii
    return D;
}

// Solves Ax=b using Jacobi's iterative method
// Matrix A must be non-singular, square and b same size as A
// It is suposed that the matrix A is already in dominant form (if possible)
// and the equation still holds. Dominant : A_{ii} > \sum_{j\neq i}{A_{ij}}
// It is necesary to provide an initial point to begin the iteration
// A good one is Vector AD = MatrixDiag(A); Vector x0 = VecDiv(b,AD);
Vector JacobiSolve(Matrix& A, Vector& b, Vector& x0, bool print, double tolerance, bool out){
    int n=A.size(); int cols=A[0].size();  // Size of M
    int bsize = b.size();                  // Size of b
    if (cols != n || bsize != n){                           // Check if M is square
        cout << "NOT SQUARE OR NOT SAME SIZE" << endl;      // error
        return {};}                                         // break
    Vector x1(n,0);                 // Definition for the next iteration
    Vector Deltax(n,0);             // Difference between iterations
    double acout[3*n+2]{};
    double epsilon = 2*tolerance;   // |x^n+1-x^n|<\epsilon
    int k{};    // initialize iteration count
    while (epsilon > tolerance){
        x1 = b;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (j != i){
                    x1[i] -= A[i][j]*x0[j];
                    x1[i] = x1[i]/A[i][i];
            }}}
        Deltax = VecDiff(x1, x0);
        epsilon = norm(Deltax);
        if (out){
            acout[0] = double(k);
            copy(x0.begin(), x0.end(), acout+1);
            copy(x1.begin(), x1.end(), acout+n+1);
            copy(Deltax.begin(), Deltax.end(), acout+n+n+1);
            acout[3*n+1] = epsilon;
            d_wa_file_Arr_csv("JacobiResult.csv", acout, 3*n+2);}
        x0 = x1;
        k++;    // advance iteration
    }
    return x1;   // Return the solution
}

Vector GaussSeidelSolve(Matrix& A, Vector& b, Vector& x0, bool print, double tolerance, bool out){
    int n=A.size(); int cols=A[0].size();  // Size of M
    int bsize = b.size();                  // Size of b
    if (cols != n || bsize != n){                           // Check if M is square
        cout << "NOT SQUARE OR NOT SAME SIZE" << endl;      // error
        return {};}                                         // break
    Vector x1(n,0);                 // Definition for the next iteration
    Vector Deltax(n,0);             // Difference between iterations
    double acout[3*n+2]{};
    double epsilon = 2*tolerance;   // |x^n+1-x^n|<\epsilon
    int k{};    // initialize iteration count
    while (epsilon > tolerance){
        x1 = b;
        for (int i =0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (j < i){
                    x1[i] -= A[i][j]*x1[j];
                    x1[i] = x1[i]/A[i][i];
                } else if (j > i){
                    x1[i] -= A[i][j]*x0[j];
                    x1[i] = x1[i]/A[i][i];
            }}}
        Deltax = VecDiff(x1, x0);
        epsilon = norm(Deltax);
        if (out){
            acout[0] = double(k);
            copy(x0.begin(), x0.end(), acout+1);
            copy(x1.begin(), x1.end(), acout+n+1);
            copy(Deltax.begin(), Deltax.end(), acout+n+n+1);
            acout[3*n+1] = epsilon;
            d_wa_file_Arr_csv("GaussSeidelResult.csv", acout, 3*n+2);}
        k++;    // advance iteration
    }
    return x1;   // Return the solution
}