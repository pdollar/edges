/*******************************************************************************
* Structured Edge Detection Toolbox      Version 3.0
* Copyright 2014 Piotr Dollar.  [pdollar-at-microsoft.com]
* Please email me if you find bugs, or have suggestions or questions!
* Licensed under the MSR-LA Full Rights License [see license.txt]
*******************************************************************************/
#include <mex.h>
#include <string.h>
#ifdef USEOMP
#include <omp.h>
#endif

typedef unsigned int uint;
template<typename T> inline T min( T x, T y ) { return x<y ? x : y; }

// run STICKY to refine superpixels (operates in place)
void sticky( uint *S, float *I, float *E, uint h, uint w, double *prm )
{
  // get additional parameters
  const uint maxIter=uint(prm[0]); uint nThreads=uint(prm[1]);
  const float sigs=float(prm[2]), sigx=float(prm[3]);
  const float sige=float(prm[4]), sigc=float(prm[5]);

  // initialize mus and ns
  uint m=0; for( uint x=0; x<w*h; x++ ) m=S[x]>m ? S[x] : m; m++;
  float *ns=new float[m](), *mus=new float[m*5]();
  for( uint x=0; x<w; x++ ) for( uint y=0; y<h; y++ ) {
    uint i=S[x*h+y]; ns[i]++; mus[i*5]+=x; mus[i*5+1]+=y;
    for(uint z=0; z<3; z++) mus[i*5+z+2]+=I[z*h*w+x*h+y];
  }
  for(uint i=0; i<m; i++) for(uint z=0; z<5; z++) mus[i*5+z]/=ns[i];

  // iterate moving boundaries
  uint changed=1, iter=0;
  while( changed && iter<maxIter ) {
    changed=0;
    #ifdef USEOMP
    nThreads = min(nThreads,uint(omp_get_max_threads()));
    #pragma omp parallel for num_threads(nThreads)
    #endif
    for( int xi=0; xi<int(w); xi++ ) for( uint y=0; y<h; y++ ) {
      // get local region, if all uniform continue
      uint x=xi, x0, x1, y0, y1, i, j, s=S[x*h+y], t, sBest, T[4];
      x0 = x>0 ? x-1 : x; x1 = x<w-1 ? x+1 : x;
      y0 = y>0 ? y-1 : y; y1 = y<h-1 ? y+1 : y;
      T[0]=S[x0*h+y]; T[1]=S[x1*h+y]; T[2]=S[x*h+y0]; T[3]=S[x*h+y1];
      if( s==T[0] && s==T[1] && s==T[2] && s==T[3] ) continue;
      // compute error of each sp label, store best
      float d, e, dBest=1e10f, ds[5], vs[5], es[4];
      vs[0]=float(x); vs[1]=float(y); sBest=T[0]+1;
      for( j=0; j<3; j++ ) vs[j+2]=I[j*h*w+x*h+y];
      es[0]=E[x0*h+y]; es[1]=E[x1*h+y]; es[2]=E[x*h+y0]; es[3]=E[x*h+y1];
      for( i=0; i<4; i++ ) {
        t=T[i]; if( t==sBest ) continue;
        e=1e10f; for(j=0; j<4; j++) if( T[j]!=t ) e=e<es[j] ? e : es[j];
        for(j=0; j<5; j++) { ds[j]=mus[t*5+j]-vs[j]; ds[j]*=ds[j]; }
        d = (ds[0]+ds[1])*sigx + (ds[2]+ds[3]+ds[4])*sigc - e*sige;
        d += sigs / ns[t]; if( d<dBest ) { dBest=d; sBest=t; }
      }
      // assign pixel to closest sp, update new and old sp
      if( sBest!=s ) {
        t=sBest; changed++; ns[s]--; ns[t]++; S[x*h+y]=t;
        for(j=0; j<5; j++) mus[s*5+j]=(mus[s*5+j]*(ns[s]+1)-vs[j])/ns[s];
        for(j=0; j<5; j++) mus[t*5+j]=(mus[t*5+j]*(ns[t]-1)+vs[j])/ns[t];
      }
    }
    if(0) mexPrintf("iter=%i changed=%f%%\n",iter,100.f*changed/h/w);
    iter++;
  }
  delete [] ns; delete [] mus; if(0) mexPrintf("--------------------\n");
}

// relabel superpixels (condense and enforce connectivity, operates in place)
void relabel( uint *S, uint h, uint w ) {
  uint *T=new uint[h*w], *map=new uint[h*w/2]();
  uint x, y, z, t, t1, t2, m=1;
  T[0]=m++; for( y=1; y<h; y++ ) T[y]=(S[y]==S[y-1]) ? T[y-1] : m++;
  for( x=1; x<w; x++ ) {
    z=x*h; T[z]=(S[z-h]==S[z]) ? T[z-h] : m++;
    for( y=1; y<h; y++ ) {
      z=y+x*h; T[z] = S[z-1]==S[z] ? T[z-1] : S[z-h]==S[z] ? T[z-h] : m++;
      if( T[z-1]!=T[z-h] && S[z-1]==S[z] && S[z-h]==S[z] ) {
        t1=T[z-1]; while(map[t1]) t1=map[t1];
        t2=T[z-h]; while(map[t2]) t2=map[t2];
        if( t1==t2 ) continue;
        if( t1<t2 ) { t=T[z]=t1; t1=T[z-h]; } else { t=T[z]=t2; t1=T[z-1]; }
        while(map[t1]) { t2=map[t1]; map[t1]=t; t1=t2; } map[t1]=t;
      }
    }
  }
  uint m1=0; for( t=1; t<m; t++ ) map[t] = map[t] ? map[map[t]] : m1++;
  for( x=0; x<w*h; x++ ) S[x]=map[T[x]];
  delete [] T; delete [] map;
}

// compute superpixel boundaries using 8-connected neighborhood
void boundaries( uint *T, uint *S, float *E, int h, int w, uint nThreads )
{
  // make 1-indexed copy of S
  for( int x=0; x<w*h; x++ ) T[x]=S[x]+1;
  // add 4-connectivity boundary greedily
  #ifdef USEOMP
  nThreads = min(nThreads,uint(omp_get_max_threads()));
  #pragma omp parallel for num_threads(nThreads)
  #endif
  for( int x=0; x<w; x++ ) for( int y=0; y<h; y++ ) {
    uint a=x*h+y, b=(x+1)*h+y, c=x*h+(y+1), s=S[a];
    if(x<w-1 && s!=S[b]) if(E[a]>E[b]) T[a]=0; else T[b]=0;
    if(y<h-1 && s!=S[c]) if(E[a]>E[c]) T[a]=0; else T[c]=0;
  }
  // add 8-connectivity boundary
  #ifdef USEOMP
  #pragma omp parallel for num_threads(nThreads)
  #endif
  for( int x=1; x<w-1; x++ ) for( int y=1; y<h-1; y++ ) {
    uint *L=T+x*h+y, i, t=L[0]; if(t==0) continue;
    const uint N[8]={L[-h],L[h],L[-1],L[1],L[-h-1],L[-h+1],L[h-1],L[h+1]};
    if( N[0] && N[1] && N[2] && N[3] ) continue;
    for( i=0; i<8; i++ ) if(N[i] && N[i]!=t) { L[0]=0; break; }
  }
  // remove excess boundary pixels
  #ifdef USEOMP
  #pragma omp parallel for num_threads(nThreads)
  #endif
  for( int x=1; x<w-1; x++ ) for( int y=1; y<h-1; y++ ) {
    uint *L=T+x*h+y, i, t=L[0]; if(t!=0) continue; t=L[0]=S[x*h+y]+1;
    const uint N[8]={L[-h],L[h],L[-1],L[1],L[-h-1],L[-h+1],L[h-1],L[h+1]};
    for( i=0; i<8; i++ ) if(N[i] && N[i]!=t) { L[0]=0; break; }
  }
}

// merge segments S that are separated by a weak boundary
void merge( uint *T, uint *S, float *E, int h, int w, float thr )
{
  // helper
  #define NEIGHBORS8(S) uint *L=S+x*h+y; int x0, y0, x1, y1; \
    x0 = (x==0) ? 0 : -h; x1 = (x==w-1) ? 0 : h; \
    y0 = (y==0) ? 0 : -1; y1 = (y==h-1) ? 0 : 1; \
    uint N[8]={L[x0],L[x1],L[y0],L[y1],L[x0+y0],L[x0+y1],L[x1+y0],L[x1+y1]};
  // compute m and min for each region
  int x; uint m=0; for( x=0; x<w*h; x++ ) m=S[x]>m ? S[x] : m; m++;
  float *es=new float[m]; for( x=0; x<int(m); x++ ) es[x]=1000;
  for( x=0; x<w*h; x++ ) es[S[x]]=min(E[x],es[S[x]]);
  // check for regions to merge and compute label mapping
  uint *map = new uint[m]();
  for( x=0; x<w; x++ ) for( int y=0; y<h; y++ ) if( S[x*h+y]==0 ) {
    uint i, j, s, s1, s2, k=0, U[8]; NEIGHBORS8(S);
    for( i=0; i<8; i++ ) if( (s=N[i])!=0 ) {
      for( j=0; j<k; j++ ) if(s==U[j]) break; if(j==k) U[k++]=s;
    }
    for( i=0; i<k-1; i++ ) for( j=i+1; j<k; j++ ) {
      s1=U[i]; while(map[s1]) s1=map[s1];
      s2=U[j]; while(map[s2]) s2=map[s2];
      if( s1!=s2 && E[x*h+y]-min(es[s1],es[s2])<thr ) {
        es[s1]=es[s2]=min(es[s1],es[s2]);
        if( s1<s2 ) { s=s1; s1=U[j]; } else { s=s2; s1=U[i]; }
        while(map[s1]) { s2=map[s1]; map[s1]=s; s1=s2; } map[s1]=s;
      }
    }
  }
  // perform mapping and remove obsolete boundaries
  uint m1=1; for( uint s=1; s<m; s++ ) map[s] = map[s] ? map[map[s]] : m1++;
  for( x=0; x<w*h; x++ ) if(S[x]) T[x]=map[S[x]];
  for( x=0; x<w; x++ ) for( int y=0; y<h; y++ ) if( T[x*h+y]==0 ) {
    uint i, s=0; NEIGHBORS8(T); for(i=0; i<8; i++) if(N[i]) { s=N[i]; break; }
    for( i; i<8; i++ ) if(N[i] && N[i]!=s) break; if(i==8) T[x*h+y]=s;
  }
  delete [] es; delete [] map;
}

// compute visualization of superpixels
void visualize( float *V, float *I, uint *S, uint h, uint w, bool bnds ) {
  uint i, z, n=h*w;
  uint m=0; for( uint x=0; x<w*h; x++ ) m=S[x]>m ? S[x] : m; m++;
  float *clrs=new float[m]; uint *cnts=new uint[n];
  for( i=0; i<m; i++ ) cnts[i]=0;
  for( i=0; i<n; i++ ) cnts[S[i]]++;
  for( z=0; z<3; z++ ) {
    for( i=0; i<m; i++ ) clrs[i]=0;
    for( i=0; i<n; i++ ) clrs[S[i]]+=I[z*n+i];
    for( i=0; i<m; i++ ) clrs[i]/=cnts[i];
    if( bnds ) clrs[0]=0;
    for( i=0; i<n; i++ ) V[z*n+i]=clrs[S[i]];
  }
  delete [] clrs; delete [] cnts;
}

// inteface to various superpixel helpers
void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] )
{
  int f; char action[1024]; f=mxGetString(pr[0],action,1024); nr--; pr++;
  if(f) { mexErrMsgTxt("Failed to get action.");

  } else if(!strcmp(action,"sticky")) {
    // S = sticky( S, I, E, prm )
    uint *S0 = (uint*) mxGetData(pr[0]);
    float *I = (float*) mxGetData(pr[1]);
    float *E = (float*) mxGetData(pr[2]);
    double *prm = (double*) mxGetData(pr[3]);
    uint h = (uint) mxGetM(pr[0]);
    uint w = (uint) mxGetN(pr[0]);
    pl[0] = mxCreateNumericMatrix(h,w,mxUINT32_CLASS,mxREAL);
    uint* S = (uint*) mxGetData(pl[0]); memcpy(S,S0,h*w*sizeof(uint));
    sticky(S,I,E,h,w,prm); relabel(S,h,w);

  } else if(!strcmp(action,"boundaries")) {
    // S = boundaries( S, E, nThreads )
    uint *S0 = (uint*) mxGetData(pr[0]);
    float *E = (float*) mxGetData(pr[1]);
    uint nThreads = (int) mxGetScalar(pr[2]);
    uint h = (uint) mxGetM(pr[0]);
    uint w = (uint) mxGetN(pr[0]);
    pl[0] = mxCreateNumericMatrix(h,w,mxUINT32_CLASS,mxREAL);
    uint* S = (uint*) mxGetData(pl[0]);
    boundaries(S,S0,E,h,w,nThreads);

  } else if(!strcmp(action,"merge")) {
    // S = wsMerge( S, E, thr );
    uint *S = (uint*) mxGetData(pr[0]);
    float *E = (float*) mxGetData(pr[1]);
    float thr = (float) mxGetScalar(pr[2]);
    uint h = (uint) mxGetM(pr[0]);
    uint w = (uint) mxGetN(pr[0]);
    pl[0] = mxCreateNumericMatrix(h,w,mxUINT32_CLASS,mxREAL);
    uint* T = (uint*) mxGetData(pl[0]);
    merge(T,S,E,h,w,thr);

  } else if(!strcmp(action,"visualize")) {
    // V = visualize( S, I, boundaries )
    uint *S = (uint*) mxGetData(pr[0]);
    float *I = (float*) mxGetData(pr[1]);
    bool bnds = mxGetScalar(pr[2])>0;
    uint h = (uint) mxGetM(pr[0]);
    uint w = (uint) mxGetN(pr[0]);
    const int dims[3] = {h,w,3};
    pl[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);
    float* V = (float*) mxGetData(pl[0]);
    visualize(V,I,S,h,w,bnds);
  }
}
