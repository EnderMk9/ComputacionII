// requires <vector> and <math.h>
// requires "rwlib.hpp" ( <fstream> and <string> )
typedef vector<vector<double>> Matrix;    // Create a type called Matrix that is a vector of vectors
typedef vector<double> Vector;            // Create a type called Vector that is a vector

//-----------------------------------------------------------------------------------------
// Cout tools
//-----------------------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------------------
// Matrix and vector operations
//-----------------------------------------------------------------------------------------

// Multiplies two matrices
Matrix matprod(Matrix& A, Matrix& B){  // Inputs are Matrices A and B passed with & as a pointer
    int Arows=A.size(); int Acols=A[0].size();  // sizes of A
    int Brows=B.size(); int Bcols=B[0].size();  // sizes of B
    if (Acols != Brows){                        // compatible multiplication check
        cout << "WRONG MULTIPLICATION" << endl; // error
        return Matrix {0};                      // break
        }
    int Crows=Arows; int Ccols=Bcols;           // dimensions of the product
    Matrix C( Crows,vector<double>(Ccols,0));   // declaration of the product
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

//-----------------------------------------------------------------------------------------
// Creation and conversion tools
//-----------------------------------------------------------------------------------------

// Creates an n-vector with all components with a certain value
Vector VecFull(double value, int n){
    Vector v(n,value);
    return v;
}

// Creates a matrix of n rows and m columns with all components having the same value
Matrix MatFull(double value, int n, int m){
    Matrix M(n,vector<double>(m,value));
    return M;
}

Matrix MIdentity(int n){
    Matrix I(n,vector<double>(n,0));
    for (int i = 0; i < n; i++){
        I[i][i] = 1;
    }
    return I;
}

Matrix Vec2Mat(Vector& V, bool col = 0){
    int n = V.size();
    Matrix M(1,vector<double>(n,0));
    M[0] = V;
    if (col){
        M = transpose(M);
    }
    return M;
}

// Creates a diagonal matrix with a vector
Matrix MDiag(Vector& V){
    int n = V.size();
    Matrix D(n,vector<double>(n,0));
    for (int i = 0; i < n; i++){
        D[i][i] = V[i];
    }
    return D;
}

// Returns the diagonal of a matrix as a vector
// A must be square
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

//-----------------------------------------------------------------------------------------
// Metric tools
//-----------------------------------------------------------------------------------------

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

// Calculates the bilinear form of two vectors defined by A
double Bilinear(Vector& X,Vector& Y, Matrix& A){
    int nX = X.size(); int nY = Y.size();
    int rows=A.size(); int cols=A[0].size();        // Size of M
    if (cols != rows || nX != rows || nY != cols){  // Check if M is square
        cout << "NOT RIGHT SIZE" << endl;           // error
        return 0;                                   // break
        }
    Matrix x = Vec2Mat(X, 0);  // X row matrix 
    Matrix y = Vec2Mat(Y, 1);  // Y column matriz
    Matrix a = matprod(x,A);   // xA
    Matrix p = matprod(a,y);   // (xA)y
    return p[0][0];
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

//-----------------------------------------------------------------------------------------
// Lineal system of equations solve, Ax = b
//-----------------------------------------------------------------------------------------

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

// a,b,c are the diagonals of the matrix as vectors
// Vector b must be 1 dimension bigger than a and c
// a and c must have the same dimension
// b and f must have the same dimension
Matrix TrDiag(Vector& a,Vector& b,Vector& c){
    int Sa=a.size(); int Sb=b.size(); int Sc=c.size(); // sizes of vectors
    if (Sa != Sc || Sb - 1 != Sa){ // Check if the sizes are correct
        cout << "INCORRECT SIZE" << endl;      // error
        return {};}  
    Matrix D( Sb,vector<double>(Sb,0));
    for (int i = 0; i < Sb-1; i++){
        D[i][i]   = b[i];
        D[i][i+1] = c[i];
        D[i+1][i] = a[i];
    }
    D[Sb-1][Sb-1] = b[Sb-1];
    return D;
}

// Solves Ax=f using Tridiagonal decomposition
// a,b,c are the diagonals of the matrix as vectors
// Vector b must be 1 dimension bigger than a and c
// a and c must have the same dimension
// b and f must have the same dimension
Vector TrDiagSolve(Vector& a,Vector& b,Vector& c, Vector& f, bool out = 0, string wname = "TROut.dat"){
    int Sa=a.size(); int Sb=b.size(); int Sc=c.size(); int Sf=f.size(); // sizes of vectors
    if (Sa != Sc || Sb - 1 != Sa || Sb != Sf){ // Check if the sizes are correct
        cout << "INCORRECT SIZE" << endl;      // error
        return {};}                                         // break
    // coutmat(L); coutmat(U);
    Vector alpha(Sa,0); Vector beta(Sb,0);// Define the vectors for the decomposition
    beta[0] = b[0];
    for (int i = 1; i < Sb; i++){
        alpha[i-1]=a[i-1]/beta[i-1];
        beta[i] = b[i]-alpha[i-1]*c[i-1];
    }
    if (out){        
            double arrbeta[Sb]{}; 
            double arralpha[Sb]{}; arralpha[Sb-1]=NAN;
            copy(alpha.begin(), alpha.end(), arralpha);
            copy(beta.begin(), beta.end(), arrbeta);
            d_w_file_2cols(wname, arralpha,arrbeta,Sb);
        }
    //coutvec(alpha); coutvec(beta);
    Vector z(Sf,0); Vector x(Sf,0);       // Define the vectors to solve
    z[0]=f[0];
    for(int i = 1; i < Sf; i++){         // First solve Lz=f
        z[i] = f[i] - alpha[i-1]*z[i-1];
      }
    //coutvec(z);
    x[Sf-1] = z[Sf-1]/beta[Sb-1];
    for(int i = Sb-2; i >= 0; i--){      // Solve Ux=z
        x[i] = (z[i]-c[i]*x[i+1])/beta[i];
        }
    return x;   // Return the solution
}

// Solves Ax=b using Jacobi's or Gauss-Seidel's iterative method
// Matrix A must be non-singular, square and b same size as A
// It is suposed that the matrix A is already in dominant form (if possible)
// and the equation still holds. Dominant : A_{ii} > \sum_{j\neq i}{A_{ij}}
// It is necesary to provide an initial point to begin the iteration
// A good one is Vector AD = MatrixDiag(A); Vector x0 = VecDiv(b,AD);
Vector JacobiSolve(Matrix& A, Vector& b, Vector& x0, double tolerance, bool print = 0){
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
        //cout << "k: " << k << endl;
        x1 = b;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (j != i){
                    x1[i] -= A[i][j]*x0[j];
                }}
            //cout << "i: " << i << endl;
            x1[i] = x1[i]/A[i][i];
        }
        Deltax = VecDiff(x1, x0);
        epsilon = norm(Deltax);
        //cout << epsilon << endl;
        if (print){
            acout[0] = double(k);
            copy(x0.begin(), x0.end(), acout+1);
            copy(x1.begin(), x1.end(), acout+n+1);
            copy(Deltax.begin(), Deltax.end(), acout+n+n+1);
            acout[3*n+1] = epsilon;
            d_wa_file_csv("JacobiResult.csv", acout, 3*n+2);}
        x0 = x1;
        k++;    // advance iteration
    }
    //cout << "Success!" << endl;
    //coutvec(x1);
    return x1;   // Return the solution
}

Vector GaussSeidelSolve(Matrix& A, Vector& b, Vector& x0, double tolerance, bool print = 0){
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
                } else if (j > i){
                    x1[i] -= A[i][j]*x0[j];
            }}
            x1[i] = x1[i]/A[i][i];
        }
        Deltax = VecDiff(x1, x0);
        epsilon = norm(Deltax);
        if (print){
            acout[0] = double(k);
            copy(x0.begin(), x0.end(), acout+1);
            copy(x1.begin(), x1.end(), acout+n+1);
            copy(Deltax.begin(), Deltax.end(), acout+n+n+1);
            acout[3*n+1] = epsilon;
            d_wa_file_csv("GaussSeidelResult.csv", acout, 3*n+2);}
        k++;    // advance iteration
    }
    return x1;   // Return the solution
}
