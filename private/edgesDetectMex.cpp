/*******************************************************************************
* Structured Edge Detection Toolbox      Version 3.01
* Code written by Piotr Dollar, 2014.
* Licensed under the MSR-LA Full Rights License [see license.txt]
*******************************************************************************/
#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef USEOMP
#include <omp.h>
#endif

typedef unsigned int uint32;
typedef unsigned short uint16;
typedef unsigned char uint8;
template<typename T> inline T min( T x, T y ) { return x<y ? x : y; }

// construct lookup array for mapping fids to channel indices
uint32* buildLookup( int *dims, int w ) {
  int c, r, z, n=w*w*dims[2]; uint32 *cids=new uint32[n]; n=0;
  for(z=0; z<dims[2]; z++) for(c=0; c<w; c++) for(r=0; r<w; r++)
    cids[n++] = z*dims[0]*dims[1] + c*dims[0] + r;
  return cids;
}

// construct lookup arrays for mapping fids for self-similarity channel
void buildLookupSs( uint32 *&cids1, uint32 *&cids2, int *dims, int w, int m ) {
  int i, j, z, z1, c, r; int locs[1024];
  int m2=m*m, n=m2*(m2-1)/2*dims[2], s=int(w/m/2.0+.5);
  cids1 = new uint32[n]; cids2 = new uint32[n]; n=0;
  for(i=0; i<m; i++) locs[i]=uint32((i+1)*(w+2*s-1)/(m+1.0)-s+.5);
  for(z=0; z<dims[2]; z++) for(i=0; i<m2; i++) for(j=i+1; j<m2; j++) {
    z1=z*dims[0]*dims[1]; n++;
    r=i%m; c=(i-r)/m; cids1[n-1]= z1 + locs[c]*dims[0] + locs[r];
    r=j%m; c=(j-r)/m; cids2[n-1]= z1 + locs[c]*dims[0] + locs[r];
  }
}

// [E,ind,segs] = mexFunction(model,I,chns,chnsSs) - helper for edgesDetect.m
void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] )
{
  // get inputs
  mxArray *model = (mxArray*) pr[0];
  float *I = (float*) mxGetData(pr[1]);
  float *chns = (float*) mxGetData(pr[2]);
  float *chnsSs = (float*) mxGetData(pr[3]);

  // extract relevant fields from model and options
  float *thrs = (float*) mxGetData(mxGetField(model,0,"thrs"));
  uint32 *fids = (uint32*) mxGetData(mxGetField(model,0,"fids"));
  uint32 *child = (uint32*) mxGetData(mxGetField(model,0,"child"));
  uint8 *segs = (uint8*) mxGetData(mxGetField(pr[0],0,"segs"));
  uint8 *nSegs = (uint8*) mxGetData(mxGetField(pr[0],0,"nSegs"));
  uint16 *eBins = (uint16*) mxGetData(mxGetField(model,0,"eBins"));
  uint32 *eBnds = (uint32*) mxGetData(mxGetField(model,0,"eBnds"));
  mxArray *opts = mxGetField(model,0,"opts");
  const int shrink = (int) mxGetScalar(mxGetField(opts,0,"shrink"));
  const int imWidth = (int) mxGetScalar(mxGetField(opts,0,"imWidth"));
  const int gtWidth = (int) mxGetScalar(mxGetField(opts,0,"gtWidth"));
  const int nChns = (int) mxGetScalar(mxGetField(opts,0,"nChns"));
  const int nCells = (int) mxGetScalar(mxGetField(opts,0,"nCells"));
  const uint32 nChnFtrs = (uint32) mxGetScalar(mxGetField(opts,0,"nChnFtrs"));
  const int stride = (int) mxGetScalar(mxGetField(opts,0,"stride"));
  const int nTreesEval = (int) mxGetScalar(mxGetField(opts,0,"nTreesEval"));
  int sharpen = (int) mxGetScalar(mxGetField(opts,0,"sharpen"));
  int nThreads = (int) mxGetScalar(mxGetField(opts,0,"nThreads"));
  const int nBnds = int(mxGetNumberOfElements(mxGetField(model,0,"eBnds"))-1)/
    int(mxGetNumberOfElements(mxGetField(model,0,"thrs")));
  const char *msgSharpen="Model supports sharpening of at most %i pixels!\n";
  if( sharpen>nBnds-1 ) { sharpen=nBnds-1; mexPrintf(msgSharpen,sharpen); }

  // get dimensions and constants
  const mwSize *imgSize = mxGetDimensions(pr[1]);
  const int h = (int) imgSize[0];
  const int w = (int) imgSize[1];
  const int Z = mxGetNumberOfDimensions(pr[1])<=2 ? 1 : imgSize[2];
  const mwSize *fidsSize = mxGetDimensions(mxGetField(model,0,"fids"));
  const int nTreeNodes = (int) fidsSize[0];
  const int nTrees = (int) fidsSize[1];
  const int h1 = (int) ceil(double(h-imWidth)/stride);
  const int w1 = (int) ceil(double(w-imWidth)/stride);
  const int h2 = h1*stride+gtWidth;
  const int w2 = w1*stride+gtWidth;
  const int imgDims[3] = {h,w,Z};
  const int chnDims[3] = {h/shrink,w/shrink,nChns};
  const int indDims[3] = {h1,w1,nTreesEval};
  const int outDims[3] = {h2,w2,1};
  const int segDims[5] = {gtWidth,gtWidth,h1,w1,nTreesEval};

  // construct lookup tables
  uint32 *iids, *eids, *cids, *cids1, *cids2;
  iids = buildLookup( (int*)imgDims, gtWidth );
  eids = buildLookup( (int*)outDims, gtWidth );
  cids = buildLookup( (int*)chnDims, imWidth/shrink );
  buildLookupSs( cids1, cids2, (int*)chnDims, imWidth/shrink, nCells );

  // create outputs
  pl[0] = mxCreateNumericArray(3,outDims,mxSINGLE_CLASS,mxREAL);
  float *E = (float*) mxGetData(pl[0]);
  pl[1] = mxCreateNumericArray(3,indDims,mxUINT32_CLASS,mxREAL);
  uint32 *ind = (uint32*) mxGetData(pl[1]);
  if(nl>2) pl[2] = mxCreateNumericArray(5,segDims,mxUINT8_CLASS,mxREAL);
  uint8 *segsOut; if(nl>2) segsOut = (uint8*) mxGetData(pl[2]);

  // apply forest to all patches and store leaf inds
  #ifdef USEOMP
  nThreads = min(nThreads,omp_get_max_threads());
  #pragma omp parallel for num_threads(nThreads)
  #endif
  for( int c=0; c<w1; c++ ) for( int t=0; t<nTreesEval; t++ ) {
    for( int r0=0; r0<2; r0++ ) for( int r=r0; r<h1; r+=2 ) {
      int o = (r*stride/shrink) + (c*stride/shrink)*h/shrink;
      // select tree to evaluate
      int t1 = ((r+c)%2*nTreesEval+t)%nTrees; uint32 k = t1*nTreeNodes;
      while( child[k] ) {
        // compute feature (either channel or self-similarity feature)
        uint32 f = fids[k]; float ftr;
        if( f<nChnFtrs ) ftr = chns[cids[f]+o]; else
          ftr = chnsSs[cids1[f-nChnFtrs]+o]-chnsSs[cids2[f-nChnFtrs]+o];
        // compare ftr to threshold and move left or right accordingly
        if( ftr < thrs[k] ) k = child[k]-1; else k = child[k];
        k += t1*nTreeNodes;
      }
      // store leaf index and update edge maps
      ind[ r + c*h1 + t*h1*w1 ] = k;
    }
  }

  // compute edge maps (avoiding collisions from parallel executions)
  if( !sharpen ) for( int c0=0; c0<gtWidth/stride; c0++ ) {
    #ifdef USEOMP
    #pragma omp parallel for num_threads(nThreads)
    #endif
    for( int c=c0; c<w1; c+=gtWidth/stride ) {
      for( int r=0; r<h1; r++ ) for( int t=0; t<nTreesEval; t++ ) {
        uint32 k = ind[ r + c*h1 + t*h1*w1 ];
        float *E1 = E + (r*stride) + (c*stride)*h2;
        int b0=eBnds[k*nBnds], b1=eBnds[k*nBnds+1]; if(b0==b1) continue;
        for( int b=b0; b<b1; b++ ) E1[eids[eBins[b]]]++;
        if(nl>2) memcpy(segsOut+(r+c*h1+t*h1*w1)*gtWidth*gtWidth,
          segs+k*gtWidth*gtWidth,gtWidth*gtWidth*sizeof(uint8));
      }
    }
  }

  // computed sharpened edge maps, snapping to local color values
  if( sharpen ) {
    // compute neighbors array
    const int g=gtWidth; uint16 N[4096*4];
    for( int c=0; c<g; c++ ) for( int r=0; r<g; r++ ) {
      int i=c*g+r; uint16 *N1=N+i*4;
      N1[0] = c>0 ? i-g : i; N1[1] = c<g-1 ? i+g : i;
      N1[2] = r>0 ? i-1 : i; N1[3] = r<g-1 ? i+1 : i;
    }
    #ifdef USEOMP
    #pragma omp parallel for num_threads(nThreads)
    #endif
    for( int c=0; c<w1; c++ ) for( int r=0; r<h1; r++ ) {
      for( int t=0; t<nTreesEval; t++ ) {
        // get current segment and copy into S
        uint32 k = ind[ r + c*h1 + t*h1*w1 ];
        int m = nSegs[k]; if( m==1 ) continue;
        uint8 S0[4096], *S=(nl<=2) ? S0 : segsOut+(r+c*h1+t*h1*w1)*g*g;
        memcpy(S,segs+k*g*g, g*g*sizeof(uint8));
        // compute color model for each segment using every other pixel
        int ci, ri, s, z; float ns[100], mus[1000];
        const float *I1 = I+(c*stride+(imWidth-g)/2)*h+r*stride+(imWidth-g)/2;
        for( s=0; s<m; s++ ) { ns[s]=0; for( z=0; z<Z; z++ ) mus[s*Z+z]=0; }
        for( ci=0; ci<g; ci+=2 ) for( ri=0; ri<g; ri+=2 ) {
          s = S[ci*g+ri]; ns[s]++;
          for( z=0; z<Z; z++ ) mus[s*Z+z]+=I1[z*h*w+ci*h+ri];
        }
        for(s=0; s<m; s++) for( z=0; z<Z; z++ ) mus[s*Z+z]/=ns[s];
        // update segment S according to local color values
        int b0=eBnds[k*nBnds], b1=eBnds[k*nBnds+sharpen];
        for( int b=b0; b<b1; b++ ) {
          float vs[10], d, e, eBest=1e10f; int i, sBest=-1, ss[4];
          for( i=0; i<4; i++ ) ss[i]=S[N[eBins[b]*4+i]];
          for( z=0; z<Z; z++ ) vs[z]=I1[iids[eBins[b]]+z*h*w];
          for( i=0; i<4; i++ ) {
            s=ss[i]; if(s==sBest) continue;
            e=0; for( z=0; z<Z; z++ ) { d=mus[s*Z+z]-vs[z]; e+=d*d; }
            if( e<eBest ) { eBest=e; sBest=s; }
          }
          S[eBins[b]]=sBest;
        }
        // convert mask to edge maps (examining expanded set of pixels)
        float *E1 = E + c*stride*h2 + r*stride; b1=eBnds[k*nBnds+sharpen+1];
        for( int b=b0; b<b1; b++ ) {
          int i=eBins[b]; uint8 s=S[i]; uint16 *N1=N+i*4;
          if( s!=S[N1[0]] || s!=S[N1[1]] || s!=S[N1[2]] || s!=S[N1[3]] )
            E1[eids[i]]++;
        }
      }
    }
  }

  // free memory
  delete [] iids; delete [] eids;
  delete [] cids; delete [] cids1; delete [] cids2;
}
