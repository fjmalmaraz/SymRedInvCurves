#define SYMMETRIC 0
#define HALF_PERIODIC 1
#define INICIO 0
#define CAMBIO 1
#define CONSERVA 2
#include <math.h>
#include <complex.h>
#include <fftw3.h>
typedef struct {
  int mod;
  int type; /* type=0 are syymetric type=1 are 1/2-periodic */
  complex *coef;
  double *tab;
  double *tabom;
  } fourier;

int dimensiona(int flag,fourier *sol, int dmod){
 sol->mod=dmod;
 if(flag==CAMBIO) { free(sol->coef);free(sol->tab); free(sol->tabom);}
 if ((sol->coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(1);}
 if ((sol->tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(1);}
 if ((sol->tabom = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(1);}   
 return 0;
}


int tabulacion( fourier *sol,double *tab){
 int j,dmod=(sol->mod);
 fftw_complex *in, *out;
 fftw_plan p;
 if((in = fftw_malloc(sizeof(fftw_complex) * dmod))==NULL){printf("no memory in tabulacion\n"); exit(1);}
 if((out = fftw_malloc(sizeof(fftw_complex) * dmod))==NULL){printf("no memory in tabulacion\n"); exit(1);}
 if(sol->type==SYMMETRIC)
  if(tab==(sol->tabom)) for(j=0;j<dmod/2;j++){ 
     in[j]=cexp(-2.*M_PI*(2.*j+1.)*I*OMEGA)*(sol->coef)[j]; in[dmod-j-1]=conj(in[j]);
  } 
  else for(j=0;j<dmod/2;j++){ 
     in[j]=(sol->coef)[j];
     in[dmod-j-1]=conj(in[j]);
  }
 else
  if(tab==(sol->tabom)) {for (j=1;j<dmod/2;j++) {in[j]=(sol->coef[j])*cexp(-4.*M_PI*j*OMEGA*I);
  in[dmod-j]= conj(in[j]);}
  in[0]=(sol->coef)[0];in[dmod/2]=0.;}
  else  {
     for(j=1;j<dmod/2 ;j++) {in[j]=(sol->coef)[j];in[dmod-j]=conj(in[j]);}
     in[dmod/2]=0.;in[0]=(sol->coef)[0];
  }
 p = fftw_plan_dft_1d(dmod, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
 fftw_execute(p);
 if (sol->type==SYMMETRIC) for(j=0;j<dmod;j++) tab[j]=creal(cexp(-M_PI*j*I/dmod)*out[j]);
 else for(j=0;j<dmod;j++) {tab[j]=creal(out[j]);} 
 fftw_destroy_plan(p);
 fftw_free(in);fftw_free(out);
 return 0;
}

int coeficiente(fourier *sol){
 int j,dmod=(sol->mod);
 fftw_complex *in, *out;
 fftw_plan p;
 if((in = fftw_malloc(sizeof(fftw_complex) * dmod))==NULL){printf("no memory in f\n"); exit(1);}
 if((out = fftw_malloc(sizeof(fftw_complex) * dmod))==NULL){printf("no memory in f\n"); exit(1);}
 if(sol->type==SYMMETRIC) for (j=0;j<dmod;j++) in[j]=cexp(M_PI*j*I/dmod)*(sol->tab)[j];
 else for(j=0;j<dmod;j++) in[j]=(sol->tab)[j]+0.*I;
 p = fftw_plan_dft_1d(dmod, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
 fftw_execute(p);
 for(j=0;j<dmod/2;j++) (sol->coef)[j]=out[j]/dmod;
 fftw_destroy_plan(p);
 fftw_free(in);fftw_free(out);
 return 0;
}

int coef_integral(fourier *sol){
 int j,dmod=(sol->mod);
 fftw_complex *in, *out;
 fftw_plan p;
 if((in = fftw_malloc(sizeof(fftw_complex) * dmod))==NULL){printf("no memory in f\n"); exit(1);}
 if((out = fftw_malloc(sizeof(fftw_complex) * dmod))==NULL){printf("no memory in f\n"); exit(1);}
 for (j=0;j<dmod;j++) in[j]=cexp(-M_PI*j*I/dmod)*(sol->tab)[j];
 p = fftw_plan_dft_1d(dmod, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
 fftw_execute(p);
 for(j=0;j<dmod/2;j++) (sol->coef)[j]=out[j]/dmod;
 fftw_destroy_plan(p);
 fftw_free(in);fftw_free(out);
 return 0;
}
