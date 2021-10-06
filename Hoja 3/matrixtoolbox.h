// requieres <vector> and <math.h>
typedef vector< vector<double> > Matrix;
typedef vector<double> Vector;

void coutmat(Matrix& M){
    cout << endl;
    int rows=M.size(); int cols=M[0].size();
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            cout << M[i][j] << " ";       // Show M matrix
        }cout << endl;
    } cout << endl;
}

Matrix matprod(Matrix& A, Matrix& B){
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

Matrix transpose(Matrix M){
    int rows=M.size(); int cols=M[0].size();
    if (cols != rows){
        cout << "NOT SQUARE" << endl;
        return Matrix {0};
        }
    for (int i = 0; i < rows; i++){
        for (int j = 0; j <= i; j++){
            swap(M[i][j],M[j][i]);
        }
    }
    return M;
}

double dotprod(Vector& X,Vector& Y){
    int Xsize=X.size(); int Ysize=Y.size(); // Size of each vector
    if (Xsize != Ysize){                    // Of sizes are not equal
        cout << "NOT SAME SIZE" << endl;    // Error
        return 0;                           // Break
        }
    double prod = 0;                        // Definition of the output
    for (int i = 0; i < Xsize; i++){        // For every component
        prod += X[i]*Y[i];                  // We multiply the entries and sum
    }
    return prod;
}

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