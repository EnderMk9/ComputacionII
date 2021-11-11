//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

typedef vector< vector<int> > MatrixInt;

void coutmatInt(MatrixInt& M){
    cout << endl;
    int rows=M.size(); int cols=M[0].size();
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            cout << M[i][j] << " ";       // Show M matrix
        }cout << endl;
    } cout << endl;
}

MatrixInt matprodInt(MatrixInt& A, MatrixInt& B){
    int Arows=A.size(); int Acols=A[0].size();
    int Brows=B.size(); int Bcols=B[0].size();
    if (Acols != Brows){
        cout << "WRONG MULTIPLICATION" << endl;
        return MatrixInt {0};
        }
    int Crows=Arows; int Ccols=Bcols;
    MatrixInt C( Crows,vector<int>(Ccols,0));
    for (int i = 0; i < Crows; i++){
        for (int j = 0; j < Ccols; j++){
            C[i][j] = 0;
            for (int k = 0; k < Crows; k++){
                C[i][j] += A[i][k]*B[k][j];     // Matrix multiplication
            }
        }
    }
    return C;
}

MatrixInt transposeInt(MatrixInt& M){
    int rows=M.size(); int cols=M[0].size();
    if (cols != rows){
        cout << "NOT SQUARE" << endl;
        return MatrixInt {0};
        }
    for (int i = 0; i < rows; i++){
        for (int j = 0; j <= i; j++){
            swap(M[i][j],M[j][i]);
        }
    }
    return M;
}