#define SYMMETRIC 0
#define HALF_PERIODIC 1
#define INICIO 0
#define CAMBIO 1
#define CONSERVA 2



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "parameter.h"

typedef struct {
  int mod;
  int type; /* type=0 are syymetric type=1 are 1/2-periodic */
  complex *coef;
  double *tab;
  double *tabom;
  } fourier;

/* int dimensiona(int flag,fourier *sol, int dmod) */

int dimensiona(int ,fourier *, int);

/* int tabulacion( fourier *sol,double *tab) */

int tabulacion( fourier *,double *);


/* int coeficiente(fourier *sol) */

int coeficiente(fourier *);

/* int coef_integral(fourier *sol) */

int coef_integral(fourier *);
