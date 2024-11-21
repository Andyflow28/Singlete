#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Qromb.h"
#include "Qromb.c"
#include "Utils.h"
#include "Utils.c"
/*
 * 
 */
#ifndef PI
#define PI   3.14159265358979323846
#endif
float A[1000][1000],B[1000][1000];
double Func1(double);
double Func2(double);
int n;
double DOS,f,g;
double Delta0 = 24.0*sqrt(2.0);

/*****************************************/
int main(int argc, char** argv) {
   
   FILE *fp;
   int i;
   int j

    ;
  
   fp = fopen("g0015delta0=1.dat","r");
   for(i=1 ; i<=1000 ; i++)
      								 for (j=1 ; j<=1; j++){
        fscanf(fp, "%f", &A[i][1]);            // Leemos un float y lo guardamos en A en la posicion i,j
        fscanf(fp, "%f", &B[i][1]);            // Leemos un float y lo guardamos en B en la posicion i,j 
       }

   fclose(fp);
   for(n=1 ; n<=1000 ; n++)
     										  //for(m=1 ; m<=2 ; m++)
       {
       										    /*C1[n]=A[n][1]									
										   C2[n]=B[n][1];*/
           printf("%f\t    %f\n",A[n][1],B[n][1]);
           
   
       }
 										  /*for(n=1 ; n<=200 ; n++)
     										  printf("%f\n",C[n]);*/
       
    		FILE *fq;
		fq = fopen("DOS2c0g0015.dat","w");       
          
                        for(n=1 ; n<=200 ; n++){
								
                        f= qromb(Func1,0.0, PI/2.0);
			g= qromb(Func2,0.0, PI/2.0);
			DOS = (f + g)/(PI/2);
                     
                fprintf(fq,"%f \t %f\n", A[n][1], DOS);    
	    }
	    fclose(fq);


   return(0);


    return (EXIT_SUCCESS);
}
/*****************************************/
double Func1(double theta)

{

double  result, aux, aux1, rho_theta, a_theta, b, positive_theta;

	a_theta = pow(A[n][1],2.0) - pow(B[n][1],2.0) - pow(Delta0*cos(2.0*theta),2.0);
	
	b = 2*A[n][1]*B[n][1];

	rho_theta = sqrt( pow(a_theta,2.0) + pow(b,2.0) );
	
		
                        positive_theta = 1 + a_theta/rho_theta;
			
			aux = sqrt(2.0*rho_theta);
			
			aux1 = A[n][1]/aux;	
	
			result = aux1*sqrt(positive_theta);

		return result;
}

/*****************************************/
double Func2(double theta)
	{	
double  result, aux2, aux3, rho_theta, a_theta, b, negative_theta;

	a_theta = pow(A[n][1],2.0) - pow(B[n][1],2.0) - pow(Delta0*cos(2.0*theta),2.0);
	
	b = 2*A[n][1]*B[n][1];

	rho_theta = sqrt( pow(a_theta,2.0) + pow(b,2.0) );

	
			negative_theta = 1 - a_theta/rho_theta;
			
			aux2 = sqrt(2.0*rho_theta);
			
			aux3 = B[n][1]/aux2;	
	
			result = aux3*sqrt(negative_theta);
        
		return result;
}
/*****************************************/