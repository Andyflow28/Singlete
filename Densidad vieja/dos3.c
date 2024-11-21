#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "Qromb.h"
#include "Qromb.c"
#include "Utils.h"
#include "Utils.c"

#define N 100
#define EULER 0.577215 
#define PI 3.141592

#define DOUBLE 1	/* use 1 for float precision and 0 for float */
#if DOUBLE
#define float double
#endif

#define EPS 1.0e-6
#define JMAX 22
#define JMAXP (JMAX+1)
#define K 5


// Interaction constants, v1-pseudo zeros, v2-horizontal nodes
float v1=0.40;
float v2=1.00;

// tight binding parameters for the  gamma band
float Hop_t = -0.40; // nearest neighbor hopping (in meV)
float Hop_t_prime = -0.12; // next nearest neighbor hopping (in meV)
float ChemPotential = -0.4; // E-E_F (in meV)
float E; // normalised chemical potential 
float t_prime; // normalised next-neighbor hopping
float UpperLimitDx; // upper limit for BZ integration

float DoS_F; // DoS at FS
float J_over_vF (float, float); // functions
float Fun1 (float,float); //  of momentum
float Fun2 (float,float); // at the Fermi
float Omega_n;
float yy1(float x);
float yy2(float x);
float quad2d(float (*func)(float,float),float x1,float x2);


float qromo(float (*func)(float x), float a, float b,float (*choose)(float(*func)(float x),float a,float b,int n));
void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
float midpnt(float (*func)(float x), float a, float b, int n);
float qromo2(float (*func)(float x), float a, float b,float (*choose)(float(*func)(float x),float a,float b,int n));
void polint2(float xa[], float ya[], int n, float x, float *y, float *dy);
float midpnt2(float (*func)(float x), float a, float b, int n);

void nrerror(char error_text[]);
float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);

static float xsav;
static float (*nrfunc)(float,float);

int main(void)
{
		FILE *fq1, *fq2;
		int n;
		float om_min, om_max, delta_om, DOS_SC, LTB_INV, LTU_INV, IMAG, MOD;
		E = ChemPotential/Hop_t;
		t_prime = Hop_t_prime/Hop_t;
		UpperLimitDx = acos(0.5*(sqrt(1-E*t_prime)-1)/t_prime);
		DoS_F = 1/PI*quad2d(J_over_vF, 0.0, UpperLimitDx);
		om_min = 0.00;
		om_max = 2.00;
		Omega_n = om_min;
		delta_om = (om_max - om_min)/(N - 1);
		fq1 = fopen("LT_BORN_INV.dat","w");
		fq2 = fopen("LT_UNITARY_INV.dat","w");

		for (n = 1; n < N; n++) {
		      Omega_n = om_min + n*delta_om;
	              DOS_SC =1/PI*quad2d(Fun1,0.0,UpperLimitDx)/DoS_F;
		      IMAG = 1/PI*quad2d(Fun2,0.0,UpperLimitDx)/DoS_F;
                      MOD = DOS_SC*DOS_SC + IMAG*IMAG;
		             if ( MOD == 0.0){
                                LTU_INV = 1.00;
			     }
                             else{ 
	                     LTB_INV = DOS_SC;
		             LTU_INV = DOS_SC/MOD;
			     } 
		      fprintf(fq1,"%f \t %f\n", Omega_n, LTB_INV);
		      fprintf(fq2,"%f \t %f\n", Omega_n,  LTU_INV);
		      }
		fclose(fq1);
		fclose(fq2);
		return 0;
}	
	

float J_over_vF(float x, float y)

{
double aux1, z, aux2, vx, vz, cos_x, sin_x, jacobian;
        double vel, cos_z, sin_z, der_zx, factor;
 
        cos_x = cos(x);
        sin_x = sqrt(1-cos_x*cos_x);
        aux1 = ChemPotential + 2*Hop_t*cos_x;
        aux2 = Hop_t + 2*Hop_t_prime*cos_x;
        z = PI-acos(0.5*aux1/aux2);
        cos_z = cos(z);
        sin_z = sqrt(1-cos_z*cos_z);
        vx = -2*(Hop_t + 2*Hop_t_prime*cos_z)*sin_x;
        vz = -2*(Hop_t + 2*Hop_t_prime*cos_x)*sin_z;
        vel = sqrt(vx*vx + vz*vz);
        der_zx = pow(-vx/vz,2);
        factor = sqrt(1 + der_zx);
        jacobian = factor/vel;
        return jacobian;
}



float Fun1(float x,float y)

// (Jacobian/v_F)*(omega/sqrt(a_k)

{	
	float result, cos_x, cos_xh, cos_z, cos_zh, sin_x, sin_xh;
	float sin_z, sin_zh, cos_yh, gap1, gap2, factor, jacobian;
	float aux1, aux2, z, vx, vz, vel, der_zx, auxgap, a_k;

        cos_x = cos(x);
        sin_x = sqrt(1-cos_x*cos_x);
        aux1 = ChemPotential + 2*Hop_t*cos_x;
        aux2 = Hop_t + 2*Hop_t_prime*cos_x;
        z = PI-acos(0.5*aux1/aux2);
        cos_z = cos(z);
        sin_z = sqrt(1-cos_z*cos_z);
        vx = -2*(Hop_t + 2*Hop_t_prime*cos_z)*sin_x;
        vz = -2*(Hop_t + 2*Hop_t_prime*cos_x)*sin_z;
        vel = sqrt(vx*vx + vz*vz);
        der_zx = pow(-vx/vz,2);
        factor = sqrt(1 + der_zx);
        jacobian = factor/vel;
	cos_xh = sqrt(0.5*(1+cos_x));
	sin_xh = sqrt(0.5*(1-cos_x));
	cos_zh = sqrt(0.5*(1+cos_z));
	sin_zh = sqrt(0.5*(1-cos_z));
	cos_yh = cos(y/2);
	gap1 = v1*sin_x + v2*sin_xh*cos_zh*cos_yh; 
        gap2 = v1*sin_z + v2*cos_xh*sin_zh*cos_yh;
	auxgap = gap1*gap1 + gap2*gap2;
	a_k = Omega_n*Omega_n - auxgap;
	
		if ( a_k <= 0)   result = 0.00;

		else    result = jacobian*Omega_n*sqrt(1/a_k);

	return result;
}


float Fun2(float x,float y)

// (Jacobian/v_F)*(-omega/sqrt(a_k)

{
	float result, cos_x, cos_xh, cos_z, cos_zh, sin_x, sin_xh;
	float sin_z, sin_zh, cos_yh, gap1, gap2, factor, jacobian;
	float aux1, aux2, z, vx, vz, vel, der_zx, auxgap, a_k;

        cos_x = cos(x);
        sin_x = sqrt(1-cos_x*cos_x);
        aux1 = ChemPotential + 2*Hop_t*cos_x;
        aux2 = Hop_t + 2*Hop_t_prime*cos_x;
        z = PI-acos(0.5*aux1/aux2);
        cos_z = cos(z);
        sin_z = sqrt(1-cos_z*cos_z);
        vx = -2*(Hop_t + 2*Hop_t_prime*cos_z)*sin_x;
        vz = -2*(Hop_t + 2*Hop_t_prime*cos_x)*sin_z;
        vel = sqrt(vx*vx + vz*vz);
        der_zx = pow(-vx/vz,2);
        factor = sqrt(1 + der_zx);
        jacobian = factor/vel;
	cos_xh = sqrt(0.5*(1+cos_x));
	sin_xh = sqrt(0.5*(1-cos_x));
	cos_zh = sqrt(0.5*(1+cos_z));
	sin_zh = sqrt(0.5*(1-cos_z));
	cos_yh = cos(y/2);
	gap1 = v1*sin_x + v2*sin_xh*cos_zh*cos_yh; 
        gap2 = v1*sin_z + v2*cos_xh*sin_zh*cos_yh;
	auxgap = gap1*gap1 + gap2*gap2;
	a_k =  auxgap - Omega_n*Omega_n;	
	
		if ( a_k <= 0.00)  result = 0;

		else 	result = -jacobian*Omega_n*sqrt(1/a_k);

	return result;
}

//**************************
// NUMERICAL RECIPES ROUTINES FOR INTEGRATION
//**************************


float yy1(float x)
{
	float c;
	
	c=0.0;
	return c;
}
float yy2(float x)
{
	float c;
	
	c=PI;
	return c;
}

float quad2d(float (*func)(float,float),float x1,float x2)
{
        float qromo2(float (*func)(float), float a, float b,float (*choose)(float(*func)(float),float a,float b,int n));
	float midpnt2(float (*func)(float), float a, float b, int n);

	float f1(float x);

	nrfunc=func;
        return qromo2(f1,x1,x2,midpnt2);
}

float f1(float x)
{
        float qromo(float (*func)(float), float a, float b,float (*choose)(float(*func)(float),float a,float b,int n));
	float midpnt(float (*func)(float), float a, float b, int n);
        float f2(float y);
        float yy1(float),yy2(float);
 
        xsav=x;
        return qromo(f2,yy1(x),yy2(x),midpnt);
}

float f2(float y)
{
	float f;
	
	f=(*nrfunc)(xsav,y);
	return f;
}


float qromo(float (*func)(float x), float a, float b,float (*choose)(float(*func)(float x),float a,float b,int n))
 
{
  void polint(float xa[],float ya[], int n,float x,float *y, float *dy);
  void nrerror(char error_text[]);
        float ss,dss;
        float s[JMAXP],h[JMAXP+1];
        int j;
 
        h[1]=1.0;
        for (j=1;j<=JMAX;j++) {
                s[j]=(*choose)(func,a,b,j);
                if (j >= K) {
                        polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
                        if (fabs(dss) <= EPS*fabs(ss)) {
                        /* Invoke the next line to print the number of iterations.  
*/
                        /* printf("interations = %d\n",j); */
                                         return ss;}
                }
                h[j+1]=h[j]/9.0;
        }
        nrerror("Too many steps in routine qromo");
        return 0.0;
}

float qromo2(float (*func)(float x), float a, float b,float (*choose)(float(*func)(float x),float a,float b,int n))
 
{
  void polint2(float xa[],float ya[], int n,float x,float *y, float *dy);
  void nrerror(char error_text[]);
        float ss,dss;
        float s[JMAXP],h[JMAXP+1];
        int j;
 
        h[1]=1.0;
        for (j=1;j<=JMAX;j++) {
                s[j]=(*choose)(func,a,b,j);
                if (j >= K) {
                        polint2(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
                        if (fabs(dss) <= EPS*fabs(ss)) {
                        /* Invoke the next line to print the number of iterations.  
*/
                        /* printf("interations = %d\n",j); */
                                         return ss;}
                }
                h[j+1]=h[j]/9.0;
        }
        nrerror("Too many steps in routine qromo2");
        return 0.0;
}


#undef EPS
#undef JMAX
#undef JMAXP
#undef K

void polint(float xa[], float ya[], int n, float x, float *y, float *dy)
{
	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}

void polint2(float xa[], float ya[], int n, float x, float *y, float *dy)
{
	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}

/* note #undef's at end of file */
#define FUNC(x) ((*func)(x))

float midpnt(float (*func)(float x), float a, float b, int n)
{
        float x,tnm,sum,del,ddel;
        static float s;
        int it,j;
 
 
        if (n == 1) {
                return (s=(b-a)FUNC(0.5(a+b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it *= 3;
                tnm=it;
                del=(b-a)/(3.0*tnm);
                ddel=del+del;
                x=a+0.5*del;
                sum=0.0;
                for (j=1;j<=it;j++) {
                  sum += FUNC(x);
                  x +=ddel;
                  sum += FUNC(x);
                  x +=del;
                }
                return (s=(s+(b-a)*sum/tnm)/3.0);
        }
}

float midpnt2(float (*func)(float x), float a, float b, int n)
{
        float x,tnm,sum,del,ddel;
        static float s;
        int it,j;
 
 
        if (n == 1) {
                return (s=(b-a)FUNC(0.5(a+b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it *= 3;
                tnm=it;
                del=(b-a)/(3.0*tnm);
                ddel=del+del;
                x=a+0.5*del;
                sum=0.0;
                for (j=1;j<=it;j++) {
                  sum += FUNC(x);
                  x +=ddel;
                  sum += FUNC(x);
                  x +=del;
                }
                return (s=(s+(b-a)*sum/tnm)/3.0);
        }
}

#undef FUNC

#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}


void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}