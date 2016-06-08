#ifndef OMEGA
#define OMEGA .618033988749894848204586834365
#endif

#ifndef EPS_FLOQUET
#define EPS_FLOQUET 1.e-7
#endif


#ifndef EPS
#define EPS 1.e-7
#endif

#ifndef LE
#define LE 0.0
#endif

#ifndef INITIAL_A
#define INITIAL_A exp(LE)
#endif

#ifndef INITIAL_B
//#define INITIAL_B 1.9
#define INITIAL_B 0.0
#endif

#ifndef INITIAL_CONTINUATION_STEP
#define INITIAL_CONTINUATION_STEP 0.01
#endif

#ifndef BD_FILE 
#define BD_FILE "bd.dat"
#endif 

#ifndef TABU_FILE
#define TABU_FILE "tabu.dat"
#endif 


#ifndef INITIAL_SOLUTION
#define INITIAL_SOLUTION 
#define INITIAL_SOLUTION_0 0.0+0.0*I
#define INITIAL_SOLUTION_1 0.0+0.0*I
#define INITIAL_SOLUTION_2 0.0+0.0*I
#endif

//#define MODOS_MAX 3000000 // In my laptop 
#define MODOS_MAX 1000 // In my macbook 
//#define MODOS_MAX 7000000 // In my desktop computer
//#define MODOS_MAX 55000000 // In kronos
//#define MODOS_MAX 80000000

