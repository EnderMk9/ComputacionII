//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

// Requires <fstream> and <string>
#include <fstream>

// Get the number of lines of a file
int file_len(string rname){
    FILE *rfile;
    rfile = fopen(rname.c_str(), "r");
    int ln = 1; char c;
    for (c = getc(rfile); c != EOF; c = getc(rfile)){
        if (c == '\n'){ // Increment count if this character is newline
            ln++;
        }
    }
    fclose(rfile);
    return ln;
}

//------------------------------------------
// Read
//------------------------------------------

// Read a list of doubles separated by each line onto an array
Vector read_lines_double(string rname){
    ifstream rfile (rname);
    int n = file_len(rname);
    Vector rdata(n);
    double line; int i = 0;
    while (rfile >> line){
        rdata[i] = line;
        i++;
    }
    rfile.close();
    return rdata;
}

// Read a list of strings separated by each line onto an array
SVector read_lines_string(string rname){
    ifstream rfile (rname);
    int n = file_len(rname);
    SVector rdata(n);
    string line; int i = 0;
    while (rfile >> line){
        rdata[i] = line;
        i++;
    }
    rfile.close();
    return rdata;
}

Matrix read_matrix_double(string rname,int m){
    ifstream rfile (rname);
    int n = file_len(rname);
    Matrix rdata(n,Vector(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
        rfile >> rdata[i][j];
        }
    }
    rfile.close();
    return rdata;
}

//------------------------------------------
// Write
//------------------------------------------

void write_col_double(string wname, Vector& wdata, int pres = 16, bool csv = 0){
    ofstream wfile (wname);  // set read file
    wfile.precision(pres);
    int n=wdata.size();
    if (not csv){
        for (int i = 0; i < n; i++){
            wfile << wdata[i] << endl;
        }   
    }else if (csv){
        for (int i = 0; i < n; i++){
            wfile << wdata[i] << "," << endl;
        } 
    }
    wfile.close();
}

void write_mat_double(string wname, Matrix& wdata, int pres = 16, bool csv = 0){
    ofstream wfile (wname);  // set read file
    wfile.precision(pres);
    int rows=wdata.size();
    int cols=wdata[0].size();
    if (not csv){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                wfile << wdata[i][j] << "  ";
            } wfile << endl;
        }
    }else if (csv){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                wfile << wdata[i][j] << ", ";
            } wfile << endl;
        }
    }
    wfile.close();
}

// writes a double array as a line without overwriting the rest of the file
void cwrite_row_double(string wname, Vector& wdata, int pres = 16, bool csv = 0){
    ofstream wfile;
    wfile.open(wname, ios_base::app);
    wfile.precision(pres);
    int n=wdata.size();
    if (not csv){
        for (int i = 0; i < n; i++){    // for every line
            wfile << wdata[i] << "  " ;
        } wfile << endl;
    }else if (csv){
        for (int i = 0; i < n; i++){    // for every line
            wfile << wdata[i] << ", " ;
        } wfile << endl;
    }
    wfile.close();  //close the file
}