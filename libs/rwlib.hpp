//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

// Requires <fstream> and <string>
#include <string>
#include <fstream>

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

void d_w_file_ln(string wname, double rdata[],int lines, int pres = 16){
    ofstream wfile (wname);  // set read file
    wfile.precision(pres);
    for (int i = 0; i < lines; i++){
        wfile << rdata[i] << endl;
    }
    wfile.close();
}

// writes two arrays as two columns separeted by two spaces in a file,
// it is necesary that each vector has the same length and that the length
// and number of lines are specified as an input.
void d_w_file_2cols(string wname, double wdata1[], double wdata2[],int lines, int pres = 16){
    ofstream wfile (wname);  // set read file
    wfile.precision(pres);
    for (int i = 0; i < lines; i++){    // for every line
        wfile << wdata1[i] << "  " << wdata2[i] << endl; // write onto the file the contents of each array separated
    }
    wfile.close();  //close the file
}
void ld_w_file_2cols(string wname, long double wdata1[], long double wdata2[],int lines, int pres = 35){
    ofstream wfile (wname);  // set read file
    wfile.precision(pres);
    for (int i = 0; i < lines; i++){    // for every line
        wfile << wdata1[i] << "  " << wdata2[i] << endl; // write onto the file the contents of each array separated
    }
    wfile.close();  //close the file
}

// writes a double array as a line without overwriting the rest of the file
void d_wa_file_csv(string waname, double wadata[], int n, int pres = 16){
    ofstream wafile;
    wafile.open(waname, ios_base::app);
    wafile.precision(pres);
    for (int i = 0; i < n; i++){    // for every line
        wafile << wadata[i] << ", " ; // write onto the file the contents of each array separated
    } wafile << endl;
    wafile.close();  //close the file
}