/*******************************************************************************
* Structured Edge Detection Toolbox      Version 3.01
* Code written by Piotr Dollar, 2014.
* Licensed under the MSR-LA Full Rights License [see license.txt]
*******************************************************************************/
#include <mex.h>
#include <math.h>
#ifdef USEOMP
#include <omp.h>
#endif

// return I[x,y] via bilinear interpolation
inline float interp( float *I, int h, int w, float x, float y ) {
  x = x<0 ? 0 : (x>w-1.001 ? w-1.001 : x);
  y = y<0 ? 0 : (y>h-1.001 ? h-1.001 : y);
  int x0=int(x), y0=int(y), x1=x0+1, y1=y0+1;
  float dx0=x-x0, dy0=y-y0, dx1=1-dx0, dy1=1-dy0;
  return I[x0*h+y0]*dx1*dy1 + I[x1*h+y0]*dx0*dy1 +
    I[x0*h+y1]*dx1*dy0 + I[x1*h+y1]*dx0*dy0;
}

// E = mexFunction(E,O,r,s,m,nThreads)
void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] )
{
  float *E0 = (float*) mxGetData(pr[0]);    // original edge map
  float *O = (float*) mxGetData(pr[1]);     // orientation map
  int r = (int) mxGetScalar(pr[2]);         // radius for nms supr
  int s = (int) mxGetScalar(pr[3]);         // radius for supr boundaries
  float m = (float) mxGetScalar(pr[4]);     // multiplier for conservative supr
  int nThreads = (int) mxGetScalar(pr[5]);  // number of threads for evaluation

  int h=(int) mxGetM(pr[0]), w=(int) mxGetN(pr[0]);
  pl[0] = mxCreateNumericMatrix(h,w,mxSINGLE_CLASS,mxREAL);
  float *E = (float*) mxGetData(pl[0]);

  // suppress edges where edge is stronger in orthogonal direction
  #ifdef USEOMP
  nThreads = nThreads<omp_get_max_threads() ? nThreads : omp_get_max_threads();
  #pragma omp parallel for num_threads(nThreads)
  #endif
  for( int x=0; x<w; x++ ) for( int y=0; y<h; y++ ) {
    float e=E[x*h+y]=E0[x*h+y]; if(!e) continue; e*=m;
    float coso=cos(O[x*h+y]), sino=sin(O[x*h+y]);
    for( int d=-r; d<=r; d++ ) if( d ) {
      float e0 = interp(E0,h,w,x+d*coso,y+d*sino);
      if(e < e0) { E[x*h+y]=0; break; }
    }
  }

  // suppress noisy edge estimates near boundaries
  s=s>w/2?w/2:s; s=s>h/2? h/2:s;
  for( int x=0; x<s; x++ ) for( int y=0; y<h; y++ ) {
    E[x*h+y]*=x/float(s); E[(w-1-x)*h+y]*=x/float(s); }
  for( int x=0; x<w; x++ ) for( int y=0; y<s; y++ ) {
    E[x*h+y]*=y/float(s); E[x*h+(h-1-y)]*=y/float(s); }
}
