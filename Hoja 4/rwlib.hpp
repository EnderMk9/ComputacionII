// Requires <fstream> and <string>

int file_ln(string rname){
    ifstream rfile (rname);
    int ln{}; string line;
    while (rfile >> line){
        ln++;}
    rfile.close();
    return ln;
}

void d_r_file_ln(string rname, double rdata[]){
    ifstream rfile (rname);
    double line; int i = 0;
    while (rfile >> line){
        rdata[i] = line;
        i++;
    }
    rfile.close();
}

void d_w_file_ln(string wname, double rdata[],int lines){
    ofstream wfile (wname);  // set read file
    for (int i = 0; i < lines; i++){
        wfile << rdata[i] << endl;
    }
    wfile.close();
}
