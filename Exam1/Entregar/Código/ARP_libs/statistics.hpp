double mean(Vector& X){
    int n = X.size(); double sum{};
    for (int i = 0; i < n; i++){
        sum += X[i];
    }
    return sum/n;
}