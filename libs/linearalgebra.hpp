//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

// requires <vector> and <math.h>
// requires "rwlib.hpp" ( <fstream> and <string> )

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
            for (int k = 0; k < Acols; k++){    // over the common index
                C[i][j] += A[i][k]*B[k][j];     // multiplication and sum
            }
        }
    }
    return C;
}

//-----------------------------------------------------------------------------------------
// Creation and conversion tools
//-----------------------------------------------------------------------------------------

Matrix Vec2Mat(Vector& V, bool col = 0){
    int n = V.size();
    Matrix M(1,vector<double>(n,0));
    M[0] = V;
    if (col){
        M = transpose(M);
    }
    return M;
}

Matrix MIdentity(int n){
    Matrix I(n,vector<double>(n,0));
    for (int i = 0; i < n; i++){
        I[i][i] = 1;
    }
    return I;
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

// calculates the standard p-norm of an array
double norm(Vector& X, int p = 2, bool inf = 0){
    if (inf){
        int n = X.size();
        Vector absX = X;
        for (int i = 0; i < n; i++){
            absX[i] = abs(X[i]);
        }
        return *max_element(absX.begin(), absX.end());
    }else if (!inf){
        int n = X.size(); double sum{};
        for (int i = 0; i < n; i++){
            sum += pow(X[i],p);
        }
        //cout << sum << endl;
        return pow(sum,1.0/p);
    }
    return {};
}

// calculates the angle between two vectors
double angle(Vector& X,Vector& Y){ // Inputs are two vectors X and Y, passed with & as a pointer
    double cos_a = dotprod(X,Y)/(sqrt(dotprod(X,X))*sqrt(dotprod(Y,Y)));
    // definition of the cosine of the angle
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

// Calculates the determinant making use of det M = det L det U
// and the determinant of a triangular matrix is the product of 
// the elements of the diagonal, L's diagonal is full of ones so
// det M = det U = Π_i^n U_ii.
double determinant(Matrix& M){             
    int n=M.size();                        // Size of M
    Matrix L( n,vector<double>(n,0));      // Lower matrix of 0s
    Matrix U( n,vector<double>(n,0));      // Upper matrix of 0s
    LU(M,L,U);
    return DiagonalProd(U);
}

// Solves Ax=b using LU decomposition
// Matrix A must be non-singular, square and b same size as A
Vector LUSolve(Matrix& A, Vector& b, bool out = 0){
    int n=A.size(); int cols=A[0].size();  // Size of M
    int bsize = b.size();                  // Size of b
    if (cols != n || bsize != n){                           // Check if M is square
        cout << "NOT SQUARE OR NOT SAME SIZE" << endl;      // error
        return {};}                                         // break
    Matrix L( n,vector<double>(n,0));      // Lower matrix of 0s
    Matrix U( n,vector<double>(n,0));      // Upper matrix of 0s
    LU(A,L,U);                             // LU decomposition
    if (out){
        cout << "L = "; coutmat(L);
        cout << "U = "; coutmat(U);
    }
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
    Vector alpha(Sb,0); Vector beta(Sb,0);// Define the vectors for the decomposition
    beta[0] = b[0];
    for (int i = 1; i < Sb; i++){
        alpha[i-1]=a[i-1]/beta[i-1];
        beta[i] = b[i]-alpha[i-1]*c[i-1];
    }
    if (out){        
            Matrix data = MatFull(0, 2, Sb);
            data[0] = alpha; data[1] = beta;
            data = transpose(data);
            write_mat_double("STrRes.dat", data);
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
Vector JacobiSolve(Matrix& A, Vector& b, Vector x0, double tol, int& k, int p = 2, bool inf = 0){
    int n=A.size(); int cols=A[0].size();  // Size of M
    int bsize = b.size();                  // Size of b
    if (cols != n || bsize != n){                           // Check if M is square
        cout << "NOT SQUARE OR NOT SAME SIZE" << endl;      // error
        return {};}                                         // break
    Vector x1(n,0);                 // Definition for the next iteration
    Vector Deltax(n,0);             // Difference between iterations
    double epsilon = 2*tol;   // |x^n+1-x^n|<\epsilon
    while (epsilon > tol){
        k++;    // advance iteration
        x1 = b;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (j != i){
                    x1[i] -= A[i][j]*x0[j];
                }}
            x1[i] = x1[i]/A[i][i];
        }
        Deltax = VecDiff(x1, x0);
        epsilon = norm(Deltax,p,inf);
        x0 = x1;
    }
    return x1;   // Return the solution
}

Vector GaussSeidelSolve(Matrix& A, Vector& b, Vector x0, double tol, int& k,int p = 2, bool inf = 0){
    int n=A.size(); int cols=A[0].size();  // Size of M
    int bsize = b.size();                  // Size of b
    if (cols != n || bsize != n){                           // Check if M is square
        cout << "NOT SQUARE OR NOT SAME SIZE" << endl;      // error
        return {};}                                         // break
    Vector x1(n,0);                 // Definition for the next iteration
    Vector Deltax(n,0);             // Difference between iterations
    double epsilon = 2*tol;   // |x^n+1-x^n|<\epsilon
    while (epsilon > tol){
        k++;    // advance iteration
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
        epsilon = norm(Deltax,p,inf);
        x0 = x1;
    }
    return x1;   // Return the solution
}

//-----------------------------------------------------------------------------------------
// Espectral theory
//-----------------------------------------------------------------------------------------

double JacobiAngle(double& aii, double& ajj, double& aij){
    if (aii != ajj){
        return atan(2*aij/(aii-ajj))/2;
    }else if(aii == ajj){
        return M_PI/4;
    }
    return {};
}

/* Matrix Rotation(double& theta, int i, int j, int n){
    Matrix R = MIdentity(n);
    R[i][i] =  cos(theta);
    R[j][j] =  cos(theta);
    R[i][j] = -sin(theta);
    R[j][i] =  sin(theta);
    return R;
} */

Vector PartMatVec(Matrix& A, IMatrix& index, bool absv = 0){
    int n = A.size();
    int p = (n*n-n)/2;
    Vector v(p,0);
    IMatrix Index(p,IVector(2,0));
    index = Index; int k{};
    if (absv){
        for (int i = 0; i < n; i++){
            for (int j = i+1; j < n; j++){
                v[k] = abs(A[i][j]);
                index[k][0] = i; index[k][1] = j;
                k++;
            }
        }
    } else if (!absv){
        for (int i = 0; i < n; i++){
            for (int j = i+1; j < n; j++){
                v[k] = A[i][j];
                index[k][0] = i; index[k][1] = j;
                k++;
            }
        }
    }
    return v;
}

IVector MaxMat(Matrix& A, double& max){
    IMatrix Index;
    Vector v = PartMatVec(A, Index,1);
    max = *max_element(v.begin(), v.end());
    int i = max_element(v.begin(), v.end())-v.begin();
    return Index[i];
}

Matrix JacobiDiag(Matrix A, Matrix& U, int& k, double tol){
    k = 0; double max;
    int n = A.size(); int i; int j;
    U = MIdentity(n); Matrix P = U;
    IVector MaxIn(2,0); double theta;
    Matrix B = A; MaxIn = MaxMat(A, max);
    while(max > tol){
        i = MaxIn[0]; j = MaxIn[1]; k++;
        theta = JacobiAngle(A[i][i],A[j][j],A[i][j]);
        for (int l = 0; l<n; l++){
            P[l][i] =  U[l][i]*cos(theta)+U[l][j]*sin(theta);
            P[l][j] = -U[l][i]*sin(theta)+U[l][j]*cos(theta);
        }
        U = P;
        for (int l = 0; l<n; l++){
            if (l!=i && l!= j){
                B[i][l] =  A[i][l]*cos(theta)+A[j][l]*sin(theta); B[l][i] = B[i][l];
                B[j][l] = -A[i][l]*sin(theta)+A[j][l]*cos(theta); B[l][j] = B[j][l];
            }
        }
        B[i][i]=A[i][i]*pow(cos(theta),2)+A[j][j]*pow(sin(theta),2)+2*A[i][j]*sin(theta)*cos(theta);
        B[j][j]=A[i][i]*pow(sin(theta),2)+A[j][j]*pow(cos(theta),2)-2*A[i][j]*sin(theta)*cos(theta);
        B[j][i] = 0; B[i][j] = 0;
        MaxIn = MaxMat(B, max);
        A = B;
    }
    return B;
}