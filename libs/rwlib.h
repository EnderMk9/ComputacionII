// Requires <fstream> and <string>

// Get the number of lines of a file
int file_ln(string rname){
    ifstream rfile (rname);
    int ln{}; string line;
    while (rfile >> line){
        ln++;}
    rfile.close();
    return ln;
}

// Read a list of doubles separated in each line onto an array
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

// writes two arrays as two columns separeted by two spaces in a file,
// it is necesary that each vector has the same length and that the length
// and number of lines are specified as an input.
void d_w_file_2cols(string wname, double rdata1[], double rdata2[],int lines){
    ofstream wfile (wname);  // set read file
    for (int i = 0; i < lines; i++){    // for every line
        wfile << rdata1[i] << "  " << rdata2[i] << endl; // write onto the file the contents of each array separated
    }
    wfile.close();  //close the file
}
void ld_w_file_2cols(string wname, long double rdata1[], long double rdata2[],int lines){
    ofstream wfile (wname);  // set read file
    for (int i = 0; i < lines; i++){    // for every line
        wfile << rdata1[i] << "  " << rdata2[i] << endl; // write onto the file the contents of each array separated
    }
    wfile.close();  //close the file
}