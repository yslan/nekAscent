#if !defined(nekascent_hpp_)
#define nekascent_hpp_

#include <vector>
#include <tuple>

#include <ascent.hpp>
#include <iostream>
#include <mpi.h>

#ifdef __cplusplus
extern "C"
{
#endif

/* the following macro functions like a##b,
   but will expand a and/or b if they are themselves macros */
#define TOKEN_PASTE_(a,b) a##b
#define TOKEN_PASTE(a,b) TOKEN_PASTE_(a,b)

#ifdef PREFIX
#  define PREFIXED_NAME(x) TOKEN_PASTE(PREFIX,x)
#else
#  define PREFIXED_NAME(x) x
#endif

#ifdef FPREFIX
#  define FPREFIXED_NAME(x) TOKEN_PASTE(FPREFIX,x)
#else
#  define FPREFIXED_NAME(x) x
#endif

#if defined(UPCASE)
#  define FORTRAN_NAME(low,up) FPREFIXED_NAME(up)
#  define FORTRAN_UNPREFIXED(low,up) up
#elif defined(UNDERSCORE)
#  define FORTRAN_NAME(low,up) FPREFIXED_NAME(TOKEN_PASTE(low,_))
#  define FORTRAN_UNPREFIXED(low,up) TOKEN_PASTE(low,_)
#else
#  define FORTRAN_NAME(low,up) FPREFIXED_NAME(low)
#  define FORTRAN_UNPREFIXED(low,up) low
#endif

//#define ascent_init     FORTRAN_UNPREFIXED(ascent_init,     ASCENT_INIT)
#define fnekascent_setup    FORTRAN_UNPREFIXED(fnekascent_setup,    FNEKASCENT_SETUP)
#define fnekascent_update   FORTRAN_UNPREFIXED(fnekascent_update,   FNEKASCENT_UPDATE)
#define fnekascent_finalize FORTRAN_UNPREFIXED(fnekascent_finalize, FNEKASCENT_FINALIZE)

typedef long long int hlong;
typedef std::vector< std::tuple<std::string, double*> > fields;

//void ascent_init();
void nekascent_setup(MPI_Comm comm_,
                     int ndim_, int lx1_, int ly1_, int lz1_, 
                     int nelt_, int lelt_, int nfldt_,
                     double *xm1_, double *ym1_, double *zm1_, 
                     double *vx, double *vy, double *vz, 
                     double *pm1, double *t);
void nekascent_update(double time, int istep, int movingMesh);
void nekascent_finalize();

void fnekascent_setup(int *comm_in,
                      int *ndim_, int *lx1_, int *ly1_, int *lz1_,             
                      int *nelt_, int *lelt_, int *nfldt_,
                      double *xm1_, double *ym1_, double *zm1_,
                      double *vx, double *vy, double *vz,
                      double *pm1, double *t);
void fnekascent_update(double *time, int *istep, int *movingMesh);
void fnekascent_finalize();


#ifdef __cplusplus
} // extern "C"
#endif
#endif
