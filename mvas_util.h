#ifndef MVAS_UTIL_INCLUDE
#define MVAS_UTIL_INCLUDE
#include <stdarg.h>
#include "nr3/cholesky.h"
#define RANI(x) (Int)((Doub)rand()/RAND_MAX*x)
#define SRAND(x) srand(x); rand()
#define DET(x) (x.nrows()==1)?x[0][0]:(x[1][1]*x[0][0]-x[1][0]*x[0][1])

inline Doub LOGDET(MatDoub m)
{
  Cholesky mchol(m);
  return mchol.logdet();
}

#define LOG(x) log(x)
#define SQR(x) (x)*(x)
#define SQRT(x) sqrt(x)
#define EXP(x) exp(x)
#define POW(x,y) (Int)pow((Doub)(x),(Doub)(y))

#define pcp(x) printf("CP> %d\n",x)
#define pcpst(x) printf("[%d] start>\n",x)
#define pcpen(x) printf("[%d] end# \n",x)
#define pcp_s(s,x) printf(" %s > %.4lf\n",s,(double)x)
//#define myp(templt, ...) fprintf(stderr,templt, ##__VAR_ARGS__)
#if defined(NDEBUG)
/* gcc's cpp has extensions; it allows for macros with a variable number of
   arguments. We use this extension here to preprocess pmesg away. */
#define pmsg(level, outfp, format, args...) ((void)0)
#else
#define DBG 5
#define LV4 4
#define PRG 0
#define ALL 10
#define CURLEV DBG
inline void pmsg(int level, FILE *outfp, const char* format, ...)
{
  va_list args;

  if (level > CURLEV)
    return;

  va_start(args, format);
  vfprintf(outfp, format, args);
  va_end(args);
}
/* print a message, if it is considered significant enough.
      Adapted from [K&R2], p. 174 */
#endif

#define SUCC 0
#define FAIL -1
#define FALSE 0
#define TRUE 1

#define takeSamp(a,b,d) for(int _i=0; _i<d; _i++) b[_i]=a[_i]

#define pvec(vec)\
    for(Int _j=0; _j<vec.size(); _j++)\
    {\
        pmsg(LV4,stdout,"%.4f ",(double)vec[_j]);\
    }\
    pmsg(LV4,stdout,"\n");

#define pvec_r(vec,vsize)\
    for(Int _j=0; _j<vsize; _j++)\
    {\
        pmsg(LV4,stdout,"%.4f ",(double)vec[_j]);\
    }\
    pmsg(LV4,stdout,"\n");



#define pvec_s(_s,vec)\
    pmsg(LV4,stdout,"%s:\n",_s);\
    for(Int _j=0; _j<vec.size(); _j++)\
    {\
        pmsg(LV4,stdout,"%.4f ",(double)vec[_j]);\
    }\
    pmsg(LV4,stdout,"\n");



#define pmat_r(mat,_r)\
    for(Int _i=0; _i<_r; _i++)\
    {\
        for(Int _j=0; _j<_r; _j++)\
        {\
            pmsg(LV4,stdout,"%.4f ",(double)mat[_i][_j]);\
        }\
        pmsg(LV4,stdout,"\n");\
    }\
 
#define pmat_pos(mat)\
    for(Int _i=0; _i<mat.nrows(); _i++)\
    {\
        for(Int _j=0; _j<mat.ncols(); _j++)\
        {\
            if(mat[_i][_j]<0) break;\
            pmsg(LV4,stdout,"%.4f ",(double)mat[_i][_j]);\
        }\
        pmsg(LV4,stdout,"\n");\
    }\
 


#define pmat(mat)\
    for(Int _i=0; _i<mat.nrows(); _i++)\
    {\
        for(Int _j=0; _j<mat.ncols(); _j++)\
        {\
            pmsg(LV4,stdout,"%.4f ",(double)mat[_i][_j]);\
        }\
        pmsg(LV4,stdout,"\n");\
    }\
 
#define pmat_s(_s,mat)\
    pmsg(LV4,stdout,"%s:\n",_s);\
    for(Int _i=0; _i<mat.nrows(); _i++)\
    {\
        for(Int _j=0; _j<mat.ncols(); _j++)\
        {\
            pmsg(LV4,stdout,"%.4f ",(double)mat[_i][_j]);\
        }\
        pmsg(LV4,stdout,"\n");\
    }\
 
inline Int maxi(VecDoub v)
{
  Doub mval = 0;
  Int mi = -1;
  Int vs = v.size();
  VecInt vm(vs);
  Int vcur = 0;
  for(Int i = 0; i < vs; i++) {
    if (mval < v[i]) {
      vcur = 0;
      mval = v[i];
      vm[vcur] = i;
      vcur++;
    } else if (mval == v[i]) {
      mval = v[i];
      vm[vcur] = i;
      vcur++;
    }
  }
  if (vcur > 0) {
    mi = vm[RANI(vcur)];
  }
  return mi;
}


// @ func - find the largest dot product value between two vector space.
inline Doub max_dotp(MatDoub a, MatDoub b)
{
  Doub maxv = 0;
  Int u = a.nrows();
  for(Int i = 0; i < a.ncols(); i++) {
    for(Int j = 0; j < b.ncols(); j++) {
      Doub val = 0;
      for(Int k = 0; k < u; k++) {
        val += a[k][i] * b[k][j];
      }
      if(abs(val) > maxv) {
        maxv = abs(val);
      }
    }
  }
  return maxv;
}

struct mvas_timer {
  clock_t st;
  clock_t en;
  Doub eclaps;
  mvas_timer() {
    eclaps = 0;
  }

  inline void start() {
    st = clock();
  }

  inline Doub end() {
    en = clock();
    eclaps = (Doub)(en - st) / CLOCKS_PER_SEC;
    return eclaps;
  }
};
#endif
