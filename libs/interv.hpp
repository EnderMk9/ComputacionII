Vector linspace(double x1, double x2, double n, bool dxcout = 0){
    Vector l(n,0);
    double dx = (x2-x1)/(n-1);
    if (dxcout){
        cout << "dx = " << dx << endl;
    }
    for (int i = 0; i < n; i++){
        l[i] = x1 + dx*i;
    }
    return l;
}

Vector countdiffval(Vector& a, double tol, IVector& index){
    int n = a.size();
    Vector r = VecFull(0,n);
    index[0] = 0;
    int l = 1; r[0] = a[0]; int p{};
    for (int i = 1; i < n; i++){
        Vector t = VecFull(1,l);
        for (int k = 0; k < l; k++){
            if (abs(a[i] - r[k]) > tol){
                t[k] = 1; // No hit, new
            } else if (abs(a[i] - r[k]) < tol){
                t[k] = 0; // Hit, repeated
                index[i] = k;
            }
        }
        if (VecProd(t) != 0){
            index[i] = l;
            r[l] = a[i]; l++;
        }   
    }
    Vector v = VecFull(0,l);
    for (int i = 0; i<l; i++){
        v[i] = r[i];
    }
    return v;
}