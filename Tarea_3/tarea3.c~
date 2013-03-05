#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

#define USAGE "./a.out nombre.txt"
void load_file(char*);
void guardar_datos(float*,char*, int);
void guardar_datos2(float*,char*, int);
void crear_matriz(gsl_matrix*,float*, int , int,int );
void print_matrix(gsl_matrix*,int,int);
void print_vector (gsl_vector*,int);
gsl_matrix* transpose(gsl_matrix*,int ,int);
gsl_matrix* multiply(gsl_matrix *,gsl_matrix *,int,int );
gsl_matrix* inverse(gsl_matrix*,int);
void crearArchivo(gsl_matrix*);
void crearArchivo2(gsl_matrix*);
gsl_matrix* cov_calculate(gsl_matrix*);
void find_eigens(gsl_matrix*, gsl_vector*,gsl_matrix*);
float mean(gsl_matrix*, int);
float sumatoria(gsl_matrix*, int, int, int);
gsl_matrix* cov_calculate2(gsl_matrix*);
float norma(gsl_matrix*,int);



int main(int argc, char **argv){
	
	FILE *in;
	int filas;
	int filas2;
	int columnas=3;
	int columnas2=3;
	char *filename=argv[1];
	char *filename2=argv[2];
	float *vector_datos;
	float *vector_datos2;
	float *vector_x;
	
	
	if(argc!=3){
   	 	printf("USAGE: %s\n", USAGE);
   	 	exit(1);
  	}
	
	
	load_file(filename);
	
	filas=contar_filas(filename,0);
	
	load_file(filename2);

	filas2=contar_filas(filename2,0);
	
	gsl_matrix *transpuesta = gsl_matrix_alloc (columnas, filas);
	gsl_matrix *matriz = gsl_matrix_alloc (filas, columnas);
	gsl_matrix *multiplicacion = gsl_matrix_alloc (filas, filas);
	gsl_matrix *vector = gsl_matrix_alloc (filas, 1);
	gsl_matrix *inversa = gsl_matrix_alloc (filas, filas);
	gsl_matrix *multiplicacion2 = gsl_matrix_alloc (columnas, 1);
	gsl_matrix *multiplicacion3 = gsl_matrix_alloc (columnas, 1);
	gsl_matrix *matriz2 = gsl_matrix_alloc (filas2, columnas2);
	gsl_matrix *covarianza = gsl_matrix_alloc (columnas2, columnas2);
	gsl_vector *eigenvals= gsl_vector_alloc(columnas2);
	gsl_matrix *eigenvecs= gsl_matrix_alloc(columnas2, columnas2);
	gsl_matrix *transpuesta2 = gsl_matrix_alloc (columnas2, filas2);

	vector_datos=malloc(2*filas*sizeof(float));
	
	vector_datos2=malloc(3*filas2*sizeof(float));
	
	guardar_datos(vector_datos,filename,filas);


	crear_matriz(matriz,vector_datos,filas,columnas,0);
	
	crear_matriz(vector,vector_datos,filas,1,1);
	
	transpuesta=transpose(matriz,columnas,filas);	
	
	multiplicacion=multiply(transpuesta,matriz,columnas,columnas);
	
	inversa=inverse(multiplicacion,columnas);
	
	multiplicacion2=multiply(transpuesta,vector,columnas,1);
	
	multiplicacion3=multiply(inversa,multiplicacion2,columnas,1);

	crearArchivo(multiplicacion3);

//////////////// CALCULO AUTOVECTORES////////////////////////

	guardar_datos2(vector_datos2,filename2,filas2);

	crear_matriz(matriz2,vector_datos2,filas2,columnas2,2);
	
	transpuesta2=transpose(matriz2,columnas2,filas2);
	
	covarianza = cov_calculate2(matriz2);	
	
	find_eigens(covarianza, eigenvals,eigenvecs); 
	
	crearArchivo2(eigenvecs);

	
}




//////////////////////////////////////////METODOS//////////////////////////////////////////////////

/////////////////////LOAD_FILE/////////////////////////
void load_file(char *filename){
	FILE *in;
	
  	
  	in = fopen(filename,"r"); 
  	
  	if(!in){
    	printf("problem opening the file %s\n", filename);
    	exit(1);
 	} 
 	
 	printf("el archivo obtenido es: %s \n",filename);
}

/////////////////////CONTAR_FILAS/////////////////////////

int contar_filas(char *filename,int selector)
{
	FILE *in;
	int contador=0;
	int c;
	
	if(!(in=fopen(filename, "r"))){
    printf("problem opening file %s\n", filename);
    exit(1);
 	}
 	
 	
 	
	do{
    c = fgetc(in);

	if(selector==0){
		if(c=='\n'){
	    	contador++;
	    }  
	}
	/*
	if(selector==1){
		if(c=="   "){
	    	contador++;
	    }
	}*/
	
	 }while(c!=EOF);
	
	fclose(in);
	return contador;
}

/////////////////////GUARDAR_DATOS/////////////////////////

void guardar_datos(float *vector,char *filename, int filas)
{
	FILE *in;
	int i;
	
	if(!(in=fopen(filename, "r"))){
    printf("problem opening file %s\n", filename);
    exit(1);
 	}
 	
 	for(i=0;i<2*filas;i++)
 	{
 		fscanf(in,"%f %f", &(vector[i*2]),&(vector[i*2+1]));
 	}
}

/////////////////////GUARDAR_DATOS2/////////////////////////

void guardar_datos2(float *vector,char *filename, int filas)
{
	FILE *in;
	int i;
	
	if(!(in=fopen(filename, "r"))){
    printf("problem opening file %s\n", filename);
    exit(1);
 	}
 	
 	for(i=0;i<3*filas;i++)
 	{
 		fscanf(in,"%f %f %f", &(vector[i*3]),&(vector[i*3+1]),&(vector[i*3+2]));
 	}
}



/////////////////////CREAR_MATRIZ/////////////////////////

void crear_matriz(gsl_matrix *matriz,float *vector, int filas, int columnas,int selector)
{
  
  int n_x=filas;
  int n_y=columnas;
  int pos;
  int i,j;


  if (selector==0)
  {
	  for(i=0;i<n_x;i++){
		for(j=0;j<n_y;j++){
		  	pos = j + (n_y * i) ;/*position in the array*/	
			gsl_matrix_set (matriz, i, 0, 1);
			gsl_matrix_set (matriz, i, 1, vector[2*i]);
			gsl_matrix_set (matriz, i, 2, (vector[2*i]*vector[2*i])/2);
			
		  
		}
	  }
  }
  if(selector==1)
  {
	  for(i=0;i<n_x;i++){
		
			gsl_matrix_set (matriz, i, 0, vector[2*i+1]);

	  }
  }
  if(selector==2)
  {
  	for(i=0;i<n_x;i++){
		for(j=0;j<n_y;j++){

			gsl_matrix_set (matriz, i, 0, vector[3*i]);
			gsl_matrix_set (matriz, i, 1, vector[3*i+1]);
			gsl_matrix_set (matriz, i, 2, vector[3*i+2]);
			
		}
	  }
  }


  

}




/////////////////////MULTIPLY/////////////////////////
gsl_matrix* multiply(gsl_matrix *matriz1,gsl_matrix *matriz2,int filas,int columnas){
		gsl_matrix *multiplicacion = gsl_matrix_alloc (filas, columnas);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, matriz1,matriz2, 0.0, multiplicacion);

		return multiplicacion;
}

/////////////////////INVERSE/////////////////////////
gsl_matrix* inverse(gsl_matrix *matriz,int filas){
	
    
    int s;
    gsl_matrix *inversa = gsl_matrix_alloc (filas, filas);
    gsl_permutation * p = gsl_permutation_alloc (filas);
    
 
    gsl_linalg_LU_decomp (matriz, p, &s);    
    gsl_linalg_LU_invert (matriz, p, inversa);
 
    
  
    gsl_permutation_free (p);
     
	   
		
 
	return inversa;
		
}


/////////////////////TRANSPOSE/////////////////////////

gsl_matrix* transpose(gsl_matrix *m,int columnas,int filas){
	
	
	gsl_matrix *transpuesta = gsl_matrix_alloc (columnas, filas);
	gsl_matrix_transpose_memcpy(transpuesta,m);



	return transpuesta;

}

/////////////////////PRINT_MATRIX/////////////////////////

void print_matrix(gsl_matrix *matriz,int filas,int columnas){

	int i=0;
	int j=0;
	for (i = 0; i < filas; ++i)
         for (j = 0; j < columnas; ++j)
             printf(j==columnas-1?" %6.3f\n":" %6.3f", gsl_matrix_get(matriz,i,j));
}



/////////////////////PRINT_VECTOR/////////////////////////


void print_vector (gsl_vector * v, int filas)
{
	int i;
 	for (i = 0; i < filas; i++) /* OUT OF RANGE ERROR */
         {
           printf ("v_%d = %g\n", i, gsl_vector_get (v, i));
         }
     
       gsl_vector_free (v);
}


/////////////////////CREAR ARCHIVO/////////////////////////

void crearArchivo(gsl_matrix *respuesta){
	FILE *in;
	char filename[50]="parametros_movimiento.dat";
	int i;
	float respuestas[3];
	in = fopen(filename,"w"); 
  	if(!in){
   	 printf("problem opening the file %s\n", filename);
   	 exit(1);
 	 }  
 	 
 	for(i=0;i<3;i++){
 		respuestas[i]=gsl_matrix_get(respuesta,i,0);
 		
 	}
 	fprintf(in,"%f %f %f",respuestas[0],respuestas[1],respuestas[2]);
 	fclose(in);
}

/////////////////////CREAR ARCHIVO 2/////////////////////////
void crearArchivo2(gsl_matrix *respuesta){
	FILE *in;
	char filename[50]="vectores_propios.dat";
	int i,j,pos;
	float respuestas[3],norm;
	in = fopen(filename,"w"); 
  	if(!in){
   	 printf("problem opening the file %s\n", filename);
   	 exit(1);
 	 }  
 	 
 	for(i=0;i<3;i++){
 	norm = norma(respuesta, i);
 		for(j=0;j<3;j++){
 		pos = j + (3 * i);
 		respuestas[pos]=(gsl_matrix_get(respuesta,i,j)/norm);
 		
 		}
 		
 	}
 	fprintf(in,"%f %f %f\n%f %f %f\n%f %f %f",respuestas[0],respuestas[1],respuestas[2],respuestas[3],respuestas[4],respuestas[5],respuestas[6],respuestas[7],respuestas[8]);
 	fclose(in);
}


/////////////////////EIGENVALUES MATRIX/////////////////////////

void find_eigens(gsl_matrix *m, gsl_vector *eigenvals,gsl_matrix *eigenvecs)
{
	gsl_eigen_symmv_workspace *w= gsl_eigen_symmv_alloc(m->size1);
	gsl_eigen_symmv(m, eigenvals, eigenvecs, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eigenvals, eigenvecs,GSL_EIGEN_SORT_ABS_ASC);
	gsl_matrix_free(m);
}


/////////////////////MEAN/////////////////////////

float mean(gsl_matrix *m, int fila)
{
	gsl_vector_view a;
	int i;
	float constante = 0;
	for (i=0;i<3;i++)
	{
	constante = constante+gsl_matrix_get(m,fila,i);	
	}

	constante = constante/3;
	return constante;

}

/////////////////////SUM/////////////////////////

float sumatoria(gsl_matrix *m, int i, int j, int M)
{
	int k;
	float di,dim,dj,djm,sum=0,resta1=0,resta2=0,multi=0;
	dim=mean(m, i);
	djm=mean(m, j);
	
	for (k=0; k<M;k++)
	{
	di = gsl_matrix_get(m,i,k);
	
	dj = gsl_matrix_get(m,j,k);
	
	
	resta1 = (di-dim);
	resta2 = (dj-djm);
	sum=(resta1*resta2)+sum;	

	}
	
	sum = sum/(M-1);
	return sum;

}

/////////////////////COVARIANCE MATRIX/////////////////////////

gsl_matrix* cov_calculate2(gsl_matrix *m)
{
	int i,j;
	float sum;
	gsl_matrix *covariance = gsl_matrix_alloc (m->size2,m->size2);
	
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
		{
		sum = sumatoria(m, i, j, 3);
		gsl_matrix_set(covariance, i, j, sum);
		}
	
	}
	return covariance;
	
}

/////////////////////NORMA/////////////////////////

float norma(gsl_matrix *m, int i)
{
	int j;
	float norma=0,a;
	
	gsl_matrix *vectores_normalizados = gsl_matrix_alloc (3,3);
	
	for (j=0;j<3;j++)
	{
		a=gsl_matrix_get(m,i,j);
		a=a*a;
		norma=a+norma;
	
	}
	norma=sqrt(norma);
	return norma;
	
}

