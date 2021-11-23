// Requiere la libreria vector, y definir las siguientes lineas sin comentar

// typedef vector<vector<double>> Matrix;  // Crea una Matriz, que es un vector de vectores
//typedef vector<double> Vector;            // Crea un vector, que es un vector, coleccion de numeros

// Esta libreria esta pensada para almacenar las funciones que vaya creando y para explicar el funcionamiento de algunos comandos, para tenerlo a mano por si lo necesito

// A.size(); define el tamaño de A. Como las matrices son vectores de vectores, esto diria el numero de filas de una matriz, pero A[0].size() diria el tamaño del vector almacenado en A[0]
// esto es, el primer vector fila de la matriz A
typedef vector<double> Vector; //Create a Vector that is a vector
typedef vector<vector<double> > Matrix; //Create a type called Matrix that is a vector of vectors
typedef vector<int> Vector_N; //Create a Vector that is a vector of integers


// coutMat enseña en consola una matriz
void coutMat(Matrix& M){
	int rows=M.size(), col=M[0].size(); //Calcula las filas y columnas
	for(int i=0; i<rows; i++){
		for(int j=0; j<col; j++){
			cout<<M[i][j]<<" "; // Muestra en pantalla cada componente utilizando bucles for
		}
		cout<<endl;
	}
}
// coutVec enseña en consola un vector
void coutVec(Vector& a){
	int rows=a.size(); // Calcula el tamaño del vector
	for(int i=0; i<rows; i++){
	cout<<a[i]<<" "; // Muestra en pantalla cada componente utilizando bucles for
	}
}

void coutVec_N(Vector_N& a){
	int rows=a.size(); // Calcula el tamaño del vector
	for(int i=0; i<rows; i++){
	cout<<a[i]<<" "; // Muestra en pantalla cada componente utilizando bucles for
	}
}
// Vecsum realiza la suma vectorial entre dos vectores
Vector Vecsum(Vector& a, Vector& b){
	int n_a=a.size(), n_b=b.size(); // Calcula los tamaños de los vectores
	if(n_a != n_b){
		cout<<endl<<endl<<"Los vectores no se pueden sumar porque no tienen la misma dimension"<<endl<<endl; // Comprueba que se puede realizar la suma
		return Vector {0};
	}
	Vector c(n_a,0); //Crea un vector lleno de unos para despues llenarlo con el vector suma
	for(int i=0; i<n_a; i++){
		c[i]=a[i]+b[i]; // Calcula la suma
	}
	return c;
}

// Calcula la multiplicacion del vector V por el escalar k
Vector V_esc(Vector A, double k){
	int rows_a=A.size();
	Vector r(rows_a,0);
	for(int i=0; i<rows_a; i++){
		r[i]=A[i]*k;
	}
	return r;
}

// Calcula la multiplicacion de la matriz M por el escalar k
Matrix M_esc(Matrix A, double k){
	int rows_a=A.size(), col_a=A[0].size();
	Matrix M(rows_a,Vector(col_a,1));
	for(int i=0; i<rows_a; i++){
		for(int j=0; j<col_a; j++){
			M[i][j]=A[i][j]*k;
		}
	}
	return M;
}

// Crea una matriz llena de unos
Matrix ones(int n_f, int n_c) {Matrix M(n_f,Vector(n_c,1)); return M;}

// Vecprod realiza el producto de componentes de un vector
Vector Vecprod(Vector& a, Vector& b){
	int n_a=a.size(), n_b=b.size();
	if(n_a != n_b){
		cout<<endl<<endl<<"Los vectores no se pueden multiplicar por componentes porque no tienen la misma dimension"<<endl<<endl;
		return Vector {0};
	}
	Vector c(n_a,0);
	for(int i=0; i<n_a; i++){
		c[i]=a[i]*b[i];
	}
	return c;
}

// Dotprod realiza el producto escalar entre dos vectores
double Dotprod(Vector& a, Vector& b){
	int n_a=a.size(), n_b=b.size();
	if(n_a != n_b){
		cout<<endl<<endl<<"El producto escalar no esta definido para vectores de distinta dimension"<<endl<<endl;
		double h=0;
		return h;
	}
	double P=0;
	for(int i=0; i<n_a; i++){
		P=P+a[i]*b[i];
	}
	return P;
}

// Norm devuelve la norma de un vector
double Norm_V(Vector& a){
	double P=0; 
	int n_a=a.size();
	for(int i=0; i<n_a; i++){
		P=P+a[i]*a[i];
	}
	P=sqrt(P);
	return P;
	
}

// Escprod devuelve un vector que es la multiplicacion por un escalar del otro vector.
Vector Escprod(Vector& a,double k){
	int n_a=a.size();
	Vector c(n_a,1);
	for(int i=0; i<n_a; i++){
		c[i]=a[i]*k;
	}
	return c;
}

// Realiza la suma de dos matrices
Matrix Matsum(Matrix A, Matrix B){
	int rows_a=A.size(), col_a=A[0].size(), rows_b=B.size(), col_b=B[0].size(); // Calcula el numero de filas y columnas de cada matriz
	if(rows_a!=rows_b || col_a!=col_b){ // Comprueba que tengan mismas dimensiones
		cout<<endl<<endl<<"No se puede realizar la suma porque tienen dimensiones distintas"<<endl<<endl;
	Matrix C={{0}};
	return C;
	}
	Matrix C=ones(rows_a,col_a);
	for(int i=0; i<rows_a; i++){
		for(int j=0; j<col_a; j++){
			C[i][j]=A[i][j]+B[i][j];
		}
	}
	return C;
}

// Para definir una matriz con sus dimensiones, la linea en la que se define C
// Matprod devuelve la multiplicacion matricial de las matrices A y B
Matrix Matprod(Matrix& A, Matrix& B){// Inputs las matrices A y B
	int rows_a=A.size(), col_a=A[0].size(), rows_b=B.size(), col_b=B[0].size(); // Calcula el numero de filas y columnas de cada matriz
	if(col_a!=rows_b){ // Comprueba que la multiplicacion matricial esta bien definida
		cout<<endl<<endl<<"El numero de columnas de A y el de filas de B no coinciden por lo que la matriz no esta bien definida"<<endl<<endl;
		return Matrix {0};
	}
	Matrix C(rows_a,vector<double>(col_b,0)); // Creamos la matriz C inicializada
	
    for (int i=0; i<rows_a; i++){            // every row
        for (int j=0; j<col_b; j++){        // every column
            for (int k=0; k<rows_a; k++){    // over the common index
                C[i][j] = C[i][j]+A[i][k]*B[k][j];     // multiplication and sum
            }
        }
    }
    return C;
}

// Crea una matriz llena de unos

// transpose transpone una matriz 
Matrix transpose(Matrix &A){
	int rows_a=A.size(), col_a=A[0].size();
	Matrix C(col_a,vector<double>(rows_a,0));
	for(int i=0; i<rows_a; i++){
		for(int j=0; j<col_a; j++){
			C[j][i]=A[i][j];
		}
	}
	return C;
}

// max_n devuelve el maximo entre dos valores
// Hacer tambien lo mismo pero de un vector, y que tambien devuelva el indice del maximo




// Intercambia las filas n1 y n2 de la matriz A

void switch_f(Matrix& A, int n1, int n2){
	int rows_a=A.size(), col_a=A[0].size();
	if((n1+1)>rows_a || (n2+1)>rows_a){ // Se comprueba que las filas que se quieren cambiar existen
		cout<<endl<<endl<<"Error, se intentan cambiar filas que no existen en la matriz"<<endl<<endl;
		return ;
	}
	Vector aux1=A[n1]; // Almacena las filas y las cambia
	Vector aux2=A[n2];
	A[n2]=aux1;
	A[n1]=aux2;
}

// Intercambia las columnas c1 y c2 de la matriz A
void switch_c(Matrix& A, int n1, int n2){
	Matrix C=transpose(A);
	switch_f(C,n1,n2);
	A=transpose(C);
}

// Comprueba que dos matrices son identicas
int check_M(Matrix A, Matrix B){
	int rows_a=A.size(), col_a=A[0].size(), rows_b=B.size(), col_b=B[0].size(); // Calcula el numero de filas y columnas de cada matriz
	if(rows_a!=rows_b || col_a!=col_b){
		cout<<endl<<endl<<"Las matrices no son iguales porque tienen dimensiones distintas"<<endl<<endl;
		int p=0;
		return p;
	}
	int p=1;
	for(int i=0; i<rows_a; i++){
		for(int j=0; j<col_a; j++){
			double h1=A[i][j], h2=B[i][j];
			if(h1==h2){
				p=p*1;
			}
			else{
				p=p*0;
				cout<<endl<<endl<<"Las matrices no son iguales porque A["<<i<<"]["<<j<<"] no es igual a B["<<i<<"]["<<j<<"]"<<endl<<endl;
				return p;
			}
		}
	}
	cout<<endl<<endl<<"Las matrices son iguales"<<endl<<endl;
	return p;
}

double max_S(double a, double b){
	if(a>b){
		return a;
	}
	else{
		return b;
	}
}

// Devuelve la maxima componente de un vector y su posicion. Si hay varios numeros iguales, devuelve el primero que encuentra
double max_V(Vector v,int &n_f){
	int rows_v=v.size();
	double a=v[0];
	n_f=0;
	for(int i=1; i<rows_v;i++){ // Compara todas las componentes con la primera y almacena la mayor de ellas
		if(v[i]>a){
			a=v[i];
			n_f=i;
		}
	}
	return a;
}

// Devuelve la maxima componente de una matriz y su posicion. Si hay varios numeros iguales, devuelve el primero que encuentra
double max_M(Matrix A, int &n_f, int &n_c){
	int rows_A=A.size(), col_A=A[0].size();
	n_f=0; n_c=0;
	Vector_N n(rows_A,0);
	Vector m(rows_A,0);
	for(int i=0; i<rows_A; i++){
		int n_aux=0; //////
		m[i]=max_V(A[i],n[i]);		
	}
	double max_A=max_V(m,n_f);
	n_c=n[n_f];
	return max_A;
}

// Calcula la derivada numerica central de la funcion f en el punto x0
double derivada(double(*f)(double), double x0, double err){
	/*
	Inputs: 
	
	f          = funcion de R en R que pide valores de tipo double
	x0         = Punto en el que queremos evaluar la derivada
	err        = Precision de la derivada
	
	Outputs:
	
	der        = Derivada numerica centrada
	
	
	*/
	double num=f(x0+err)-f(x0-err);
	double der=num/(2*err);
	return der;
}

// Calcula la derivada parcial de F respecto a la i+1 variable (recordemos que la primera componente de un vector es si i=0)
double der_par(double(*F)(Vector), Vector x0, double h, int i){
	/*
	Inputs: 
	
	F          = funcion de R^n en R que pide valores de tipo Vector (variables)
	x0         = Vector en el que queremos evaluar la derivada parcial
	err        = Precision de la derivada
	i          = Variable respecto de la que queremos derivar
	
	Outputs:
	
	der_par    = Derivada numerica centrada parcial 
	
	*/
	int dim_var=x0.size(); // Calculamos el numero de variables
	if((i+1)>dim_var){ // Nos aseguramos de que queremos derivar respecto de una variable existente
		cout<<endl<<endl<<"CUIDADO"<<endl<<"Se esta intentando calcular la derivada parcial respecto de una variable que no existe"<<endl<<endl<<endl<<"CUIDADO";
	}
	// Calculamos el punto x0 con un desplazamiento h en la direccion iesima, hacia delante y hacia detras (forward y backward), y evaluamos la funcion
	Vector e_i(dim_var,0); e_i[i]=1;   // Creamos un vector en la direccion iesima 
	Vector e_h=V_esc(e_i,h);           // Lo escalamos con el factor h
	Vector e_hh=V_esc(e_h,-1);         // Calculamos su inverso
	Vector x_f=x0; x_f=Vecsum(x_f,e_h);
	Vector x_b=x0; x_b=Vecsum(x_b,e_hh);
	double F_f=F(x_f), F_b=F(x_b);
	double der_par=(F_f-F_b)/(2*h);
	return der_par;
}


// Calcula la matriz jacobiana de la funcion F
Matrix Jacobian(Vector(*F)(Vector), Vector x0, double h){
	/*
	Inputs: 
	
	F          = funcion de R^n en R^m que pide valores de tipo Vector (variables)
	x0         = Vector en el que queremos evaluar el jacobiano
	err        = Precision de las derivadas
	
	Outputs:
	
	Jacobian   = Matriz Jacobiana
	
	*/	 
	
	int n=x0.size(); //      Calculo el numero de variables, esto es la dimension del espacio de salida de F
	Vector Fx=F(x0); //      Calculo la imagen de x0
	int m=Fx.size(); //      Calculo la dimension del espacio de llegada de F
	Matrix J=ones(n,m); //   Creo la Matriz Jacobiana y la lleno de unos, para luego sustituir por el valor correspondiente
	for(int j=0; j<m; j++){
      	// Aqui trabajamos con Fx[i]
		Vector e_i(n,0); e_i[j]=1;   // Creamos un vector en la direccion j-esima 
		Vector e_h=V_esc(e_i,h);           // Lo escalamos con el factor h, para evaluar F(x0+e_h)
		Vector e_hh=V_esc(e_h,-1);         // Calculamos su inverso, para evaluar F(x0-e_h)
		Vector x_f=x0; x_f=Vecsum(x_f,e_h);
		Vector x_b=x0; x_b=Vecsum(x_b,e_hh);
		Vector F_F=F(x_f), F_B=F(x_b);
		F_B=V_esc(F_B,-1);
		Vector num=Vecsum(F_F,F_B); double den=1.0/(2*h);
		Vector fr=V_esc(num,den);
		J[j]=fr; // Aqui almaceno en la fila j-esima las componentes calculadas. Sin embargo, esto no es el jacobiano ya que el jacobiano tiene por filas los gradientes de Fi
		// Por ello hay que al calcular J, trasponerlo para que quede bien
	}		
	// Tal como esta definido el programa, es la traspuesta asi que trapones J y ya queda todo bien
	J=transpose(J);
	return J;
}





// Metodos de resolucion de ecuaciones no lineales

//                  METODO DE LA BISECCION, la buena
// Converge y es lineal, no suele dar problemas


double biseccion(double a, double b, double eps, double(*funcion)(double), int &n) {
	/*
	Inputs: 
	
	a          =  Primer valor del intervalo a estudiar
	b          =  Ultimo valor del intervalo a estudiar
	eps        =  TamaÃ±o minimo para aceptar la raiz z abs(b-c)<= eps
	funcion    =  La funcion cuyos ceros vamos a estudiar
	n          =  Numero de iteraciones del metodo
	
	Outputs:
	
	c          =  Raiz encontrada en el intervalo 
	n          =  Numero de iteraciones. No es un output como tal pero esta funcion modifica directamente el valor de n en el programa principal
	*/
	n=0;
	double fa=funcion(a),fb=funcion(b); // Calculamos la funcion en los extremos del intervalo
	if( fa*fb>=0){ // Primero comprobamos que o bien no se puede aplicar Bolzano o bien tenemos ya un 0
		if(fa==0) { // Calculamos si hay un 0 en a
			cout<<endl<<endl<<"cero en a"<<endl<<endl;
		}
		else if(fb==0){// Calculamos si hay un 0 en b
			cout<<endl<<endl<<"cero en b"<<endl<<endl;
		}
		else{// En el caso de que se cumpla este if, fa*fb es estrictamente mayor que 0, no se puede aplicar Bolzano
			cout<<endl<<endl<<"Error; En este intervalo no se puede aplicar el Teorema de Bolzano pues f(a)*f(b)>0"<<endl<<endl;
		}
	}
	else {
		double c; // Es necesario definir c ahora para que al salir del bucle do-while siga estando definida
		double aux; // Lo mismo ocurre con aux, es necesario definirlo fuera del bucle y que luego se modifique dentro
		do{
			c=(a+b)/2; // Calculamos el punto c entre medias del intervalo
			n++; // Incrementamos el contador n
			double fc=funcion(c); // Calculamos el valor de la funcion en c
			if(fc==0) {
				return(c);
			}
			// En caso de que c no sea una raiz, se dan dos situaciones, o bien el 0 esta entre a y c, o esta entre b y c. La funcion calcula esto ahora precisamente, pues sabemos que 
			// estara en el intervalo en el que los signos en los extremos sean distintos
			else if (fa*fc<0) 
			{b=c;}
			else
			{a=c;}
			// Una vez calculado el nuevo intervalo, comprobamos si c es suficientemente pequeño como para valer como una raiz aproximada
			aux=max(abs(b-c),abs(a-c));		
			}while(abs(aux)>=eps) ;
		return(c);
		}
} 


//              METODO DE LA SECANTE
//   Similar al de newton pero tiene mas interes cuando la derivada numerica supone un coste computacional alto y da mas problemas.
double secante(double x0, double x1, double eps, double(*f)(double), int &n) {
		/* 
	Inputs: 
	
	x0         =  Primer valor proximo a la raiz
	x1         =  Segundo valor proximo a la raiz
	eps        =  Tolerancia que queremos para definir la proximidad entre el valor hayado y la verdadera raiz
	f          =  La funcion cuyo cero vamos a estudiar
	contador   =  Numero de iteraciones necesarias
	
	Outputs:
	
	c          =  Raiz encontrada */
	n=0;
	if(abs(f(x0))<eps){return x0;} if(f(x1)<eps){return x1;} // Devolvemos automaticamente x0 o x1 si son raices 
	// Si no lo son, necesitamos calcular el bucle
	do {
		double sec=(f(x1)-f(x0))/(x1-x0); // Calculamos la pendiente de la secante
		x0=x1; // Redefinimos el siguiente paso, primero cambiando x0
		x1=x1-f(x1)/sec;
		n++;
	} while (abs(x1-x0)>eps);
	return(x1);
}


//                  METODO DE NEWTON-RAPHSON con derivada numerica
// Puede diverger en algunas situaciones, conviene pintar la grafica, hacerse una idea de donde esta el 0 y empezar por ahi. si ves que hace cosas raras lo paras.	
double Newt_Raph_n(double x0, double eps, double(*funcion)(double), int &n) {
		/* 
	Inputs: 
	
	x0         =  Primer valor proximo a la raiz
	b          =  Ultimo valor del intervalo a estudiar
	eps        =  Tolerancia que queremos para definir la proximidad entre el valor hayado y la verdadera raiz
	funcion    =  La funcion cuyo cero vamos a estudiar
	n          =  Numero de iteraciones del metodo
	
	Outputs:
	
	c          =  Raiz encontrada 
	n          =  Numero de iteraciones. No es un output como tal pero esta funcion modifica directamente el valor de n en el programa principal
	*/
	double x=x0;
	n=0;
	do {
		double h=eps, der_x=(funcion(x+h)-funcion(x-h))/(2*h); // Calculamos la derivada centrada en xn.
		n++; // Incrementamos el contador
		if(funcion(x)==0) {return(x);}
		else if(der_x==0 && funcion(x)!=0) {
		cout<<"Error; se ha hallado punto con derivada paralela al eje y con f(x) no nulo, no se puede continuar el metodo"<<endl<<endl<<"CUIDADO"<<endl<<endl<<"CUIDADO";
		break;
	}
		else {
			x=x-funcion(x)/der_x;
		}
		} while (abs(funcion(x))>=eps);
	return(x);
}



