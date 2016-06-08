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


double coefficient_parameter_system(fourier *tint,fourier *rf,double lambda,double dmod){
  // tint is fourier series with the coefficients calculated from the integral term 
  // rf is the forcing term 
  double sum=0.;
  int j;
  if ((rf->coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){
            printf("no memory in variable rf\n"); exit(EXIT_FAILURE);
   }
   coeficiente(rf);
   
   
   for(j=dmod/2-1;-1<j;j--){ 
     if(cabs( rf->coef[j])>1.e-16){
       sum+=creal((tint->coef)[j]*(rf->coef)[j]/(cexp(-2.*M_PI*(2*j+1)*OMEGA*I)-lambda));	 
     }
   }
     
   return 2*sum; 
}

double restarting_from_old_solution(double *a,double *b,fourier *u,double h,double aold,double bold,fourier *uold, double aold1,double bold1,fourier *uold1){
  int j;
  double hold,difa=aold-aold1,difb=bold-bold1; 
  hold=sqrt(difa*difa+difb*difb);
  *a=aold+difa*h/hold;
  *b=bold+difb*h/hold;
  for(j=0;j<(uold1->mod)/2;j++) (u->coef)[j]=(uold->coef)[j]+((uold->coef)[j]-(uold1->coef)[j])*h/hold;
  for(j=(uold1->mod)/2;j<(uold->mod)/2;j++) (u->coef)[j]=(uold->coef)[j]*(1.0+h/hold);
  for(j=(uold->mod)/2;j<(u->mod)/2;j++) (u->coef)[j]=0.;
  return EXIT_SUCCESS; 
}
     
int main(){
  int ndim=101,j,iout,iter,ifin=1,dmod=ndim-1;
  double a=INITIAL_A,sum,lon=0.,h=INITIAL_CONTINUATION_STEP,hold,aold,*Fzero,liapunov,hb,ha,b,bold,bold1,error,aux,autabj,coef1b,coef1a,coef10,coef2b,coef2a,coef20,det,aold1;
  fourier u,uold,uold1,cred,daF,dbF,itF,tint,y,tpseudoa;
  FILE *bd,*tabu;
  
  u.type=uold.type=uold1.type=dbF.type=daF.type=itF.type=tint.type=y.type=tpseudoa.type=SYMMETRIC;
  cred.type=HALF_PERIODIC; 
  printf("*** Continuation symmetric invariant curves *******\n");

  if ((Fzero = (double*) malloc(ndim*sizeof(double)))==NULL){printf("no memory in f (1)\n"); exit(EXIT_FAILURE);}
  dimensiona(INICIO,&u,dmod);
  dimensiona(INICIO,&cred,dmod);
  dimensiona(INICIO,&uold,dmod);dimensiona(INICIO,&uold1,0);
  

  b=INITIAL_B;
  bold=b-h;
  bold1=bold-h;
  aold=a; 
  aold1=aold;
  hold=h;
  u.coef[0]=INITIAL_SOLUTION_0;
  u.coef[1]=INITIAL_SOLUTION_1;
  u.coef[2]=INITIAL_SOLUTION_2;
  for (j=3;j<dmod/2;j++) u.coef[j]=0.;
  uold.coef[0]=INITIAL_SOLUTION_0;
  uold.coef[1]=INITIAL_SOLUTION_1;
  uold.coef[2]=INITIAL_SOLUTION_2;
  for (j=3;j<dmod/2;j++) uold.coef[j]=0.;
  
  printf("hinicial: %f\n", h);
 
 if((bd=fopen(BD_FILE,"w"))==NULL) { printf("Error opening file bd.dat!\n"); exit(EXIT_FAILURE); }
 if((tabu=fopen(TABU_FILE,"w"))==NULL) { printf("Error opening file tabu.dat!\n"); exit(EXIT_FAILURE); }
  while(ifin>0){
   iout=1;   iter=0;
   
   while(iout>0){
     // ERROR CALCULATION IN PREVIOUS ITERATION  
     tabulacion(&u,u.tab);
     tabulacion(&u,u.tabom);
     func(&u,Fzero, a, b); 
     error=0.;
     for(j=0;j<ndim;j++){
       aux=fabs(Fzero[j]); if(error<aux) error=aux;
     }
     aux=fabs((a-aold)*(aold-aold1)+(b-bold)*(bold-bold1)-hold*h);
     if(error<aux) error=aux;
     printf("Error: %e\n",error);
     if ((iter>5)||(error>2.)){
	h*=0.5;
        if (fabs(h)<EPS*EPS) salir( "Too small step",aold,bold,&uold);
        restarting_from_old_solution(&a,&b,&u,h,aold,bold,&uold,aold1,bold1,&uold1);
        iter=0;
        printf("No convergence for h=%e\n ",h);
       }
     else if(error< EPS){iout=0;}
       else{
         //REDUCTION
         for(j=0;j<dmod;j++) cred.tab[j]=log(a*df(a*u.tab[j]));
         coeficiente(&cred);
         liapunov=exp(creal(cred.coef[0])); 
         cred.coef[0]=0.;
         for(j=1;j<dmod/2;j++) {if(cabs((cred.coef)[j])<1.e-16) cred.coef[j]=0.;
                              else cred.coef[j]/=(cexp(-4.*M_PI*j*OMEGA*I)-1.); }
         tabulacion(&cred,cred.tab);  tabulacion(&cred,cred.tabom); 
         for(j=0;j<dmod;j++){cred.tab[j]=exp(cred.tab[j]);cred.tabom[j]=exp(cred.tabom[j]);}
         sum=fabs(liapunov*cred.tabom[0]-cred.tab[0]*a*df(a*u.tab[0])); 
         for(j=1;j<dmod;j++) {
	  aux=fabs(liapunov*cred.tabom[j]-cred.tab[j]*a*df(a*u.tab[j]));
          if( aux>sum) sum=aux;
         }
         printf("Error reduccion: %e \n",sum);
         if(sum>EPS_FLOQUET) {printf("Cambio\n");
	   dmod+=2*(dmod/10); ndim=dmod+1;iter=0; 
           free(Fzero);
           if ((Fzero=(double*)malloc(ndim*sizeof(double)))==NULL){printf("no memory in f (2)"); exit(EXIT_FAILURE);}
	   if(dmod > MODOS_MAX){ salir("Excede el tamagno",aold,bold,&uold); }
            dimensiona(CAMBIO,&u,dmod); dimensiona(CAMBIO,&cred,dmod); 
	   h*=0.8;
           restarting_from_old_solution(&a,&b,&u,h,aold,bold,&uold,aold1,bold1,&uold1);           
           printf("Fin Cambio");
         } //close condition change if not reducible 
         else{
           /***************************************************
            ***  Solve system for Newton's method           ***
            ***************************************************/
           // Define size of the vectors which define system.
	   tint.mod=dbF.mod=daF.mod=itF.mod=y.mod=tpseudoa.mod=dmod;
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
             //tint.tab[j]=a*ddf(autabj)/df(autabj)*cred.tab[j];
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
            coef1b=coefficient_parameter_system(&tint,&dbF,liapunov,dmod);
            free(dbF.tab);
            if ((itF.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(EXIT_FAILURE);}
            for(j=0;j<dmod;j++) itF.tab[j]=Fzero[j]/cred.tabom[j];	    
            coef10=coefficient_parameter_system(&tint,&itF,liapunov,dmod)-Fzero[dmod];
            free(itF.tab);
            // I need the step size for the integral is the double since later we have to multiply by two
            coef1a=0;
            for(j=0;j<dmod;j++){
              autabj=a*u.tab[j];
              coef1a+=u.tab[j]*ddf(autabj)/df(autabj)/((float)dmod);
	    }
        
            if ((daF.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(EXIT_FAILURE);}
            for(j=0;j<dmod;j++)daF.tab[j]=df(a*u.tab[j])*u.tab[j]/cred.tabom[j];	    
            coef1a+=coefficient_parameter_system(&tint,&daF,liapunov,dmod);
            coef1a+=1.0/a;
            free(tint.coef);free(daF.tab);
       
        
            //if ((tpseudoa.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(EXIT_FAILURE);}

         
            //for(j=0;j<uold1.mod/2;j++) tpseudoa.tab[j]=(uold.tab[j]-uold1.tab[j])*cred.tab[j];  
            //for(j=uold1.mod/2;j<uold.mod/2;j++) tpseudoa.tab[j]=uold.tab[j]*cred.tab[j];
	    // for(j=uold.mod/2;j<u.mod/2;j++) tpseudoa.tab[j]=0.;

	    //if ((tpseudoa.coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(EXIT_FAILURE);}
            //coef_integral(&tpseudoa);
            //free(tpseudoa.tab);
        
	    //coef2b=coefficient_parameter_system(&tpseudoa,&dbF,liapunov,dmod)+bold-bold1;
	

	    //sum=0; 
            //for(j=0;j<uold1.mod;j++) sum+=(u.tab[j]-uold.tab[j])*(uold.tab[j]-uold1.tab[j]);
            //for(j=uold1.mod;j<uold.mod;j++) sum+=(u.tab[j]-uold.tab[j])*(uold.tab[j]);
            //for(j=uold.mod;j<uold.mod;j++) sum+=0;
            //sum+=;
            //sum+=(b-bold)*(bold-bold1);
            //coef20=coefficient_parameter_system(&tpseudoa,&itF,liapunov,dmod)-
            coef2a=aold-aold1;
            coef2b=bold-bold1;          
            //coef20=h*h-(a-aold)*(aold-aold1)-(b-bold)*(bold-bold1);
            coef20=hold*h-(a-aold)*(aold-aold1)-(b-bold)*(bold-bold1);
            /////CONTINUATION STEP
          
            // I need the step size for the integral is the double since later we have to multiply by two


            //coef2a=coefficient_parameter_system(&tpseudoa,&daF,liapunov,dmod)+aold-aold1;
            //SOLUTION OF THE SYSTEM CRAMER'S RULE
            printf("%f %f %f %f %f %f\n",coef1a,coef1b,coef10,coef2a,coef2b,coef20);
	    det=coef1a*coef2b-coef2a*coef1b; 
            ha=(coef10*coef2b-coef20*coef1b)/det;
            hb=(coef1a*coef20-coef2a*coef10)/det;
            printf("%f %f\n",ha,hb);
            if ((y.coef = (complex*) calloc((dmod/2),sizeof(complex)))==NULL){printf("no memory in u\n"); exit(EXIT_FAILURE);}
            for(j=0;j<dmod/2;j++) {
	      if(cabs(ha*daF.coef[j]+hb*dbF.coef[j]-itF.coef[j])>1.e-16){
	        y.coef[j]=(ha*daF.coef[j]+hb*dbF.coef[j]-itF.coef[j])/(cexp(-2.*M_PI*(2*j+1)*OMEGA*I)-liapunov);
	       }
	      else y.coef[j]=0.;
            }
            free(daF.coef);free(dbF.coef);free(itF.coef);
            if ((y.tab = (double*) calloc((dmod),sizeof(double)))==NULL){printf("no memory in tabulacion\n"); exit(EXIT_FAILURE);}
            tabulacion(&y,y.tab);
            for(j=0;j<dmod;j++) y.tab[j]*=cred.tab[j];
            coeficiente(&y);
            b+=hb; a+=ha;	iter++;
            for(j=0;j<dmod/2;j++) u.coef[j]+=y.coef[j];
            free(y.coef);free(y.tab);
	 } //Close Newton's iteration or change size for reduction (everything done if error >EPS 
       }// close iter >5
   } //Close while iout>0

   //sum=0.;
   //for(j=1;j<8;j++) {
   //       sum+=creal(u.coef[u.mod/2-j]*conj(u.coef[u.mod/2-j]));
   //}
   sum=0;
   for(j=1;j<dmod/10;j++) {
     aux=creal(u.coef[u.mod/2-j]*conj(u.coef[u.mod/2-j]));
     if(sum<aux) sum=aux;
   } 
     printf("para cambio: %e\n",sum);
   if(sum>EPS*EPS*EPS){ 
     printf("Cambio\n");dmod+=2*(dmod/10); ndim=dmod+1;iter=0;
     
     
     free(Fzero);
     if ((Fzero=(double*)malloc(ndim*sizeof(double)))==NULL){printf("no memory in f(3)"); exit(EXIT_FAILURE);}
      if(dmod > MODOS_MAX)salir("Too much memory needed",aold,bold,&uold);
     
      dimensiona(CAMBIO,&u,dmod);      
      dimensiona(CAMBIO,&cred,dmod);
     restarting_from_old_solution(&a,&b,&u,h,aold,bold,&uold,aold1,bold1,&uold1); 
     
     printf("Fin Cambio\n");
   }
   else{
     //OUTPUT SOLUTION 
     printf("*****Iteration: %d dimension: %d  a: %e error: %e\n", iter, ndim,a,error);
      sum=0.25/((double) dmod*dmod);
     lon=0.;for(j=1;j<dmod;j++) {
       aux=u.tab[j]-u.tab[j-1];
        lon += sqrt(aux*aux+sum);
     }  
     //if((a>5.)&&(h>0.2*(a-aold)/(pow(lon/lonold,3.)-1.))){ printf("Se ha passado\n"); 
     //  h=0.2*(a-aold)/(pow(lon/lonold,3.)-1.);
     //}
    
     fprintf(bd,"%e %e %e %d %e %e \n",a,b,lon,iter,h,error);
     fflush(bd); 
     tabulacion(&u,u.tab);
     for (j=0;j<dmod;j++)    fprintf(tabu,"   %e \n",u.tab[j]);
     fprintf(tabu,"\n");
     //UPDATE STEP 
     if(iter==0) h*=2.;
     if(iter==1) h*=1.5;
     if(iter==2) h*=1.1;
     if(iter==3) h*=0.75;
     if(iter>3 ) h*=0.5; 
     if (fabs(h)<1.e-09)salir("too much small step (end)",a,b,&uold);
     aux=a-aold;
     hold=aux*aux;
     aux=b-bold;
     hold+=aux*aux;
     hold=sqrt(hold);
     //UPDATING DATA FOR NEXT STEP
     aold1=aold;aold=a;
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
     restarting_from_old_solution(&a,&b,&u,h,aold,bold,&uold,aold1,bold1,&uold1); 
     
     /*
     for(j=0;j<uold1.mod/2;j++) u.coef[j]=uold.coef[j]+(uold.coef[j]-uold1.coef[j])*h/hold;
     for(j=uold1.mod/2;j<uold.mod/2;j++) u.coef[j]=uold.coef[j];
     for(j=uold.mod/2;j<u.mod/2;j++) u.coef[j]=0.; */
   } //Close if fourier series is smail queue and output solution
  } //Close while ifin>1     
  fclose(bd);
    
  return EXIT_SUCCESS;
}

