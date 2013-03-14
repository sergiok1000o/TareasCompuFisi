#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int main (int argc, char **argv)
{
	double a = atof(argv[1]);
	double b = atof(argv[2]);
	double h = atof(argv[3]);
	int i;
	int tamano;
	double k = a;
	double integral = 0.0;
	double fun, fun2;


	tamano = (int)((b - a) / h);
	

	for(i=0;i<tamano;i++)
	{
		fun = 1 + cos(k)*sin(k);
		fun2 = 1/sqrt (fun)  ;
		k = a + h;
		integral = integral + fun2 * h;
	}
 	
        printf(" el resultado de la integral es %f \n",integral);

	return 0;
}


