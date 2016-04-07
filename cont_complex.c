#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <bits/nan.h>

#include "parameter.h"
#include "cont_complex.h"

//#define MODOS_MAX 3000000 // In my laptop 
//#define MODOS_MAX 7000000 // In my desktop computer
//#define MODOS_MAX 55000000 // In kronos
#define MODOS_MAX 80000000

double f(double x){
   return atan(x);
}

double df(double x){
  return 1.0/(1.0+x*x);
}

double g(double x){
  return sin(x);
}

int func( fourier *sol,  double *fourierSeries,double a,double b){
 int j,dmod=(sol->mod);
 double sum;
 fourierSeries[dmod]=0.;
 for(j=0;j<dmod;j++){
   sum=(sol->tab)[j];
   fourierSeries[j]=((sol->tabom)[j])-f(a*sum)-b*g(M_PI*j/((double)dmod));
   fourierSeries[dmod]+=log(fabs(a*df(a*sum)));
  }
  fourierSeries[dmod]/= (double) dmod;
  fourierSeries[dmod]-=LE;
   
 return 0;
}

int salir(char *mensaje,double a,double b,fourier *u){
  FILE *sol;
  printf(mensaje);
  if((sol=fopen("sol.dat","wb"))==NULL){printf("Error opening file sol.dat!\n"); exit(1); }
  fwrite(&a,sizeof(double),1,sol);
  fwrite(&b,sizeof(double),1,sol);
  fwrite(&(u->mod),sizeof(int),1,sol);
  fwrite(u->coef,sizeof(complex),u->mod,sol);	  
  fclose(sol);
  printf(mensaje);
  exit(1);
  return 0;
  }
  
int main(){
  int ndim=101,j,iout,iter,ifin=1,dmod=ndim-1,prim;
  double a=2.44153,sum,lon=0.,h=0.01,aold=a,*fourierSeries,hold=h,liapunov,hb,b,bold,bold1,error,aux,lonold=0.;
  fourier u,uold,uold1,cred,dbF,itF,tint,y;
  FILE *bd,*tabu;
  
  u.type=uold.type=uold1.type=dbF.type=itF.type=tint.type=y.type=SYMMETRIC;
  cred.type=HALF_PERIODIC; 
  printf("Programa continua bifurcaciones\n");

  if ((fourierSeries = (double*) malloc(ndim*sizeof(double)))==NULL){printf("no memory in f (1)\n"); exit(1);}
  dimensiona(INICIO,&u,dmod);
  dimensiona(INICIO,&cred,dmod);
  dimensiona(INICIO,&uold,dmod);dimensiona(INICIO,&uold1,0);


  bold1=bold=b=1.9;
  u.coef[0]=(2.6220033115E-01-7.5384363165E-01*I)/2.;
  u.coef[1]=(-1.9204713762E-01-7.6114155828E-02*I)/2.;
  u.coef[2]=1.1766929986E-01/2.;
  for (j=3;j<dmod/2;j++) u.coef[j]=0.;
  uold.coef[0]=(2.6220033115E-01-7.5384363165E-01*I)/2.;
  uold.coef[1]=(-1.9204713762E-01-7.6114155828E-02*I)/2.;
  uold.coef[2]=1.1766929986E-01/2.;
  for (j=3;j<dmod/2;j++) uold.coef[j]=0.;
  
  
 
 if((bd=fopen("bd.dat","w"))==NULL) { printf("Error opening file bd.dat!\n"); exit(1); }
 if((tabu=fopen("tabu.dat","w"))==NULL) { printf("Error opening file bd.dat!\n"); exit(1); }
 prim=0;
  while(ifin>0){
   iout=1;   iter=0;
   while(iout>0){
     tabulacion(&u,u.tab);
     tabulacion(&u,u.tabom);
     func(&u,fourierSeries, a, b); 
     error=0.;
     for(j=0;j<ndim;j++) if(error<fabs(fourierSeries[j])) error=fabs(fourierSeries[j]);
     printf("Error: %e\n",error);
     if ((iter>5)||(error>2.)){
	h*=0.5;
        if (fabs(h)<EPS*EPS) salir( "Paso pequegn0",a,b,&u);
        a=aold+h;
        b=bold+(bold-bold1)*h/hold;
        for(j=0;j<uold1.mod/2;j++) u.coef[j]=uold.coef[j]+(uold.coef[j]-uold1.coef[j])*h/hold;
        for(j=uold1.mod/2;j<uold.mod/2;j++) u.coef[j]=uold.coef[j];
        for(j=uold.mod/2;j<u.mod/2;j++) u.coef[j]=0;
        iter=0;
        printf("No convergencia h=%e\n ",h);
       }
     else if(error< EPS){iout=0;}
     else{
       for(j=0;j<dmod;j++) cred.tab[j]=log(a*df(a*u.tab[j]));
       coeficiente(&cred);
       liapunov=exp(creal(cred.coef[0])); 
       cred.coef[0]=0.;
       for(j=1;j<dmod/2;j++) {if(cabs((cred.coef)[j])<1.e-16) cred.coef[j]=0.;
                              else cred.coef[j]/=(cexp(-4.*M_PI*j*OMEGA*I)-1.); }
       tabulacion(&cred,cred.tab);  tabulacion(&cred,cred.tabom); 
       for(j=0;j<dmod;j++){cred.tab[j]=exp(cred.tab[j]);cred.tabom[j]=exp(cred.tabom[j]);}
       sum=0.; for(j=1;j<dmod;j++) {aux=liapunov*cred.tabom[j]-cred.tab[j]*a/(1.+a*a*u.tab[j]*u.tab[j]);
                     sum+=aux*aux;}
       printf("Error reduccion: %e \n",sum);
       if(sum>EPS*EPS) {printf("Cambio\n");
	   dmod+=2*(dmod/10); ndim=dmod+1;iter=0; 
           free(fourierSeries);
           if ((fourierSeries=(double*)malloc(ndim*sizeof(double)))==NULL){printf("no memory in f (2)"); exit(1);}
	   if(dmod > MODOS_MAX){ salir("Excede el tamagno",a,b,&u); 
	   }
           dimensiona(CAMBIO,&u,dmod); ;dimensiona(CAMBIO,&cred,dmod);
	   b=bold+(bold-bold1)*h/hold;
	   for(j=0;j<uold1.mod/2;j++) u.coef[j]=uold.coef[j]+(uold.coef[j]-uold1.coef[j])*h/hold;
           for(j=uold1.mod/2;j<uold.mod/2;j++) u.coef[j]=uold.coef[j];
           for(j=uold.mod/2;j<u.mod/2;j++) u.coef[j]=0.;
           printf("Fin Cambio");
        }
       else{
         /***************************************************
          ***  Solve system for Newton's method           ***
          ***************************************************/
         // Define size of the vectors which define system.
	 tint.mod=dbF.mod=itF.mod=y.mod=dmod;
	 // tint is a vector for aprox. of integral
         /* where the vector d is the derivative of the map over the solution 
           and this is used in order to aprox. Liapunov exponent  for the next
           solution 
           L_j:= \sum_{k=0}^{n-1}  e^{-i\pi k/n} d(\theta_k) c(\theta_k) e^{-i2\pi j k /n} 
	 */ 
        if ((tint.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(1);}
        for(j=0;j<dmod;j++){
	  tint.tab[j]=-2.*a*a*u.tab[j]*cred.tab[j]*df(a*u.tab[j]);
        }
	if ((tint.coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(1);}
        coef_integral(&tint);
        free(tint.tab);
        if ((dbF.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(1);}
        for(j=0;j<dmod;j++){ 
         dbF.tab[j]=sin(M_PI*j/dmod)/cred.tabom[j]; 
        }
        if ((dbF.coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(1);}
        coeficiente(&dbF);
        free(dbF.tab);
        sum=0.; 
        for(j=dmod/2-1;-1<j;j--){ 
         if(cabs( dbF.coef[j])>1.e-16){
           sum+=creal(tint.coef[j]*dbF.coef[j]/(cexp(-2.*M_PI*(2*j+1)*OMEGA*I)-liapunov));
	 }
        }  
	sum*=2.; 
        if ((itF.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(1);}
        for(j=0;j<dmod;j++){ 
	    itF.tab[j]=fourierSeries[j]/cred.tabom[j];	    
        }
        if ((itF.coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(1);}
        coeficiente(&itF);
        free(itF.tab);
        hb=0.;
        for(j=dmod/2-1;-1<j;j--){ 
	   if(cabs( itF.coef[j])>1.e-16) hb-=tint.coef[j]*itF.coef[j]/(cexp(-2.*M_PI*(2*j+1)*OMEGA*I)-liapunov);
        }
	hb*=2.; hb+=fourierSeries[dmod]; hb/=sum;
        free(tint.coef);
        if ((y.coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(1);}
        for(j=0;j<dmod/2;j++) {if(cabs(hb*dbF.coef[j]+itF.coef[j])>1.e-16) y.coef[j]=(hb*dbF.coef[j]+itF.coef[j])/(cexp(-2.*M_PI*(2*j+1)*OMEGA*I)-liapunov);
	          else y.coef[j]=0.;}
        free(dbF.coef);free(itF.coef);
        if ((y.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(1);}
        tabulacion(&y,y.tab);
        for(j=0;j<dmod;j++) y.tab[j]*=cred.tab[j];
        coeficiente(&y);
        b-=hb; for(j=0;j<dmod/2;j++) u.coef[j]-=y.coef[j];
        free(y.coef);free(y.tab);
	iter++;

       } 
     }
   }

   sum=0.;
   for(j=1;j<8;j++) {
          sum+=creal(u.coef[u.mod/2-j]*conj(u.coef[u.mod/2-j]));
   }
   printf("para cambio: %e\n",sum);
     if(sum>EPS*EPS*EPS*8){ 
     printf("Cambio\n");dmod*=2; ndim=dmod+1;iter=0; 
     
     
     free(fourierSeries);
     if ((fourierSeries=(double*)malloc(ndim*sizeof(double)))==NULL){printf("no memory in f(3)"); exit(1);}
      if(dmod > MODOS_MAX)salir("Excede el tamagno",a,b,&u);
     
     dimensiona(CAMBIO,&u,dmod);      
     dimensiona(CAMBIO,&cred,dmod);
     
     for(j=0;j<uold1.mod/2;j++) u.coef[j]=uold.coef[j]+(uold.coef[j]-uold1.coef[j])*h/hold;
     for(j=uold1.mod/2;j<uold.mod/2;j++) u.coef[j]=uold.coef[j];
     for(j=uold.mod/2;j<u.mod/2;j++) u.coef[j]=0.;
     
     printf("Fin Cambio\n");
     }
   else{
     hold=h;
     if(iter==0) h*=2.;
     if(iter==1) h*=1.5;
     if(iter==4) h*=0.8;
     if(iter==5) h*=0.5; 
     if (fabs(h)<1.e-09)salir("too much small step (end)",a,b,&u);
    
     printf("Iteraciones: %d dimension: %d  a: %e\n", iter, ndim,a); 
     tabulacion(&u,u.tab);
     for (j=0;j<dmod;j++)    fprintf(tabu,"   %e \n",u.tab[j]);
     fprintf(tabu,"\n");
     lonold=lon;
     lon=0.;for(j=1;j<dmod;j++) lon += sqrt((u.tab[j]-u.tab[j-1])*(u.tab[j]-u.tab[j-1])+0.25/((double) dmod*dmod)); 
     if((a>5.)&&(h>0.2*(a-aold)/(pow(lon/lonold,3.)-1.))){ printf("Se ha passado\n"); 
     h=0.2*(a-aold)/(pow(lon/lonold,3.)-1.);}
    
     fprintf(bd,"%e %e %e\n",a,b,lon);
     fflush(bd); 
     aold=a;
     a+=h;
     bold1=bold; bold=b;
     if(uold1.mod<uold.mod){
       uold1.mod=uold.mod;
       free(uold1.coef);
       if ((uold1.coef = (complex*) calloc((uold1.mod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(1);}
     }
     for(j=0;j<uold1.mod/2;j++) uold1.coef[j]=uold.coef[j];
     if(uold.mod<u.mod){
       uold.mod=u.mod;
       free(uold.coef);
       if ((uold.coef = (complex*) calloc((uold.mod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(1);}
     }
     for(j=0;j<uold.mod/2;j++) uold.coef[j]=u.coef[j];
     for(j=0;j<uold1.mod/2;j++) u.coef[j]=uold.coef[j]+(uold.coef[j]-uold1.coef[j])*h/hold;
     for(j=uold1.mod/2;j<uold.mod/2;j++) u.coef[j]=uold.coef[j];
      for(j=uold.mod/2;j<u.mod/2;j++) u.coef[j]=0.;
      prim++;
    }
  }     
  fclose(bd);
    
  return 0;
}
