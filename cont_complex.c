#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>


#include "parameter.h"
#include "symmetric_curves.h"


double f(double x){
   return atan(x);
}

double df(double x){
  return 1.0/(1.0+x*x);
}

double ddf(double x){
  double aux;
  aux=df(x);
  return -2.0*aux*aux*x;
}
double g(double x){
  return sin(x);
}

int func( fourier *sol,  double *Fzero,double a,double b){
 int j,dmod=(sol->mod);
 double sum;
 Fzero[dmod]=0.;
 for(j=0;j<dmod;j++){
   sum=(sol->tab)[j];
   Fzero[j]=((sol->tabom)[j])-f(a*sum)-b*g(M_PI*j/((double)dmod));
   Fzero[dmod]+=log(fabs(a*df(a*sum)));
  }
  Fzero[dmod]/= (double) dmod;
  Fzero[dmod]-=LE;
   
 return EXIT_SUCCESS;
}

int salir(char *mensaje,double a,double b,fourier *u){
  FILE *sol;
  printf("%s\n",mensaje);
  if((sol=fopen("sol.dat","wb"))==NULL){printf("Error opening file sol.dat!\n"); exit(EXIT_FAILURE); }
  fwrite(&a,sizeof(double),1,sol);
  fwrite(&b,sizeof(double),1,sol);
  fwrite(&(u->mod),sizeof(int),1,sol);
  fwrite(u->coef,sizeof(complex),u->mod,sol);	  
  fclose(sol);
  printf("%s\n",mensaje);
  exit(EXIT_FAILURE);
  return EXIT_SUCCESS;
  }
  
int main(){
  int ndim=101,j,iout,iter,ifin=1,dmod=ndim-1,prim;
  double a=INITIAL_A,sum,lon=0.,h=0.01,aold=a,*Fzero,hold=h,liapunov,hb,b,bold,bold1,error,aux,lonold=0.,autabj;
  fourier u,uold,uold1,cred,dbF,itF,tint,y;
  FILE *bd,*tabu;
  
  u.type=uold.type=uold1.type=dbF.type=itF.type=tint.type=y.type=SYMMETRIC;
  cred.type=HALF_PERIODIC; 
  printf("Programa continua bifurcaciones\n");

  if ((Fzero = (double*) malloc(ndim*sizeof(double)))==NULL){printf("no memory in f (1)\n"); exit(EXIT_FAILURE);}
  dimensiona(INICIO,&u,dmod);
  dimensiona(INICIO,&cred,dmod);
  dimensiona(INICIO,&uold,dmod);dimensiona(INICIO,&uold1,0);


  bold1=bold=b=INITIAL_B;
  u.coef[0]=INITIAL_SOLUTION_0;
  u.coef[1]=INITIAL_SOLUTION_1;
  u.coef[2]=INITIAL_SOLUTION_2;
  for (j=3;j<dmod/2;j++) u.coef[j]=0.;
  uold.coef[0]=INITIAL_SOLUTION_0;
  uold.coef[1]=INITIAL_SOLUTION_1;
  uold.coef[2]=INITIAL_SOLUTION_2;
  for (j=3;j<dmod/2;j++) uold.coef[j]=0.;
  
  
 
 if((bd=fopen("bd.dat","w"))==NULL) { printf("Error opening file bd.dat!\n"); exit(EXIT_FAILURE); }
 if((tabu=fopen("tabu.dat","w"))==NULL) { printf("Error opening file bd.dat!\n"); exit(EXIT_FAILURE); }
 prim=0;
  while(ifin>0){
   iout=1;   iter=0;
   while(iout>0){
     tabulacion(&u,u.tab);
     tabulacion(&u,u.tabom);
     func(&u,Fzero, a, b); 
     error=0.;
     for(j=0;j<ndim;j++) if(error<fabs(Fzero[j])) error=fabs(Fzero[j]);
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
       sum=0.; for(j=1;j<dmod;j++) {aux=liapunov*cred.tabom[j]-cred.tab[j]*a*df(a*u.tab[j]);
                     sum+=aux*aux;}
       printf("Error reduccion: %e \n",sum);
       if(sum>EPS*EPS) {printf("Cambio\n");
	   dmod+=2*(dmod/10); ndim=dmod+1;iter=0; 
           free(Fzero);
           if ((Fzero=(double*)malloc(ndim*sizeof(double)))==NULL){printf("no memory in f (2)"); exit(EXIT_FAILURE);}
	   if(dmod > MODOS_MAX){ salir("Excede el tamagno",a,b,&u); }
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
        if ((tint.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(EXIT_FAILURE);}
        for(j=0;j<dmod;j++){
          autabj=a*u.tab[j];
          //For feneral functions f and g 
          //tint.tab[j]=a*ddf(autabj)/df(autabj);
          //This is the most efficient  calculation for the example in the paper 
	  tint.tab[j]=-2.*a*autabj*cred.tab[j]*df(autabj);
        }
	if ((tint.coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(EXIT_FAILURE);}
        coef_integral(&tint);
        free(tint.tab);
        if ((dbF.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(EXIT_FAILURE);}
        for(j=0;j<dmod;j++){ 
         dbF.tab[j]=g(M_PI*j/dmod)/cred.tabom[j]; 
        }
        if ((dbF.coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(EXIT_FAILURE);}
        coeficiente(&dbF);
        free(dbF.tab);
        sum=0.; 
        for(j=dmod/2-1;-1<j;j--){ 
         if(cabs( dbF.coef[j])>1.e-16){
           sum+=creal(tint.coef[j]*dbF.coef[j]/(cexp(-2.*M_PI*(2*j+1)*OMEGA*I)-liapunov));
	 }
        }  
	sum*=2.; 
        if ((itF.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(EXIT_FAILURE);}
        for(j=0;j<dmod;j++){ 
	    itF.tab[j]=Fzero[j]/cred.tabom[j];	    
        }
        if ((itF.coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(EXIT_FAILURE);}
        coeficiente(&itF);
        free(itF.tab);
        hb=0.;
        for(j=dmod/2-1;-1<j;j--){ 
	   if(cabs( itF.coef[j])>1.e-16) hb-=tint.coef[j]*itF.coef[j]/(cexp(-2.*M_PI*(2*j+1)*OMEGA*I)-liapunov);
        }
	hb*=2.; hb+=Fzero[dmod]; hb/=sum;
        free(tint.coef);
        if ((y.coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(EXIT_FAILURE);}
        for(j=0;j<dmod/2;j++) {if(cabs(hb*dbF.coef[j]+itF.coef[j])>1.e-16) y.coef[j]=(hb*dbF.coef[j]+itF.coef[j])/(cexp(-2.*M_PI*(2*j+1)*OMEGA*I)-liapunov);
	          else y.coef[j]=0.;}
        free(dbF.coef);free(itF.coef);
        if ((y.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(EXIT_FAILURE);}
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
     
     
     free(Fzero);
     if ((Fzero=(double*)malloc(ndim*sizeof(double)))==NULL){printf("no memory in f(3)"); exit(EXIT_FAILURE);}
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
     if(iter>4 ) h*=0.5; 
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
       if ((uold1.coef = (complex*) calloc((uold1.mod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(EXIT_FAILURE);}
     }
     for(j=0;j<uold1.mod/2;j++) uold1.coef[j]=uold.coef[j];
     if(uold.mod<u.mod){
       uold.mod=u.mod;
       free(uold.coef);
       if ((uold.coef = (complex*) calloc((uold.mod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(EXIT_FAILURE);}
     }
     for(j=0;j<uold.mod/2;j++) uold.coef[j]=u.coef[j];
     for(j=0;j<uold1.mod/2;j++) u.coef[j]=uold.coef[j]+(uold.coef[j]-uold1.coef[j])*h/hold;
     for(j=uold1.mod/2;j<uold.mod/2;j++) u.coef[j]=uold.coef[j];
      for(j=uold.mod/2;j<u.mod/2;j++) u.coef[j]=0.;
      prim++;
    }
  }     
  fclose(bd);
    
  return EXIT_SUCCESS;
}
