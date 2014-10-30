function [S,V] = spDetect( I, E, varargin )
% Detect Sticky Superpixels in image.
%
% Detect "Sticky Edge Adhesive Superpixels" in image. High quality, fast
% superpixels that "stick" to edges. Without edge term the code computes
% superpixels using an iterative approach motivated by both SLIC (Achanta
% et al., PAMI12) and SEEDS (Bergh et al., ECCV12) superpixels. With edge
% term added, the superpixels snap to edges, resulting in higher quality
% boundaries. There is no corresponding publication for this code at this
% time but please cite our edge detection work if you use this code.
%
% The most important parameter is k which controls superpixel scale.
% Note that the edge image E is optional (that is E=[] may be used).
%
% USAGE
%  opts = spDetect()
%  [S,V] = spDetect( I, [E], [opts] )
%
% INPUTS
%  I          - [h x w x 3] color input image (in [0,255])
%  E          - [h x w] type single edge image (in [0,1]), or [] array
%  opts       - parameters (struct or name/value pairs)
%   .type       - ['sticky'] options are 'sticky' or 'watershed'
%   .nIter      - [4] number of iterations
%   .nThreads   - [4] number of computation threads
%   .k          - [512] controls scale of superpixels (big k -> big sp)
%   .alpha      - [.5] relative importance of regularity versus data terms
%   .beta       - [.9] relative importance of edge versus color terms
%   .merge      - [0] set to small value to merge nearby superpixels at end
%   .bounds     - [1] if true add boundaries to superpixels
%   .seed       - [] optional initial seed superpixels
%
% OUTPUTS
%  S          - [h x w] superpixel label map (S==0 are boundaries)
%  V          - [h x w] superpixel visualization
%
% EXAMPLE
%
% See also spDemo, spAffinities, watershed
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

% get default parameters
dfs = { 'type','sticky', 'nIter',4, 'nThreads',4, 'k',512, ...
  'alpha',.5, 'beta',.9, 'merge',0, 'bounds',1, 'seed',[] };
o = getPrmDflt(varargin,dfs,1); if(nargin==0), S=o; return; end
type=lower(o.type(1)); assert( type=='w' || type=='s' );
sigs = [ o.k*o.alpha/1e4 o.alpha/1e4 ...
  (1-o.alpha)*o.beta (1-o.alpha)*(1-o.beta) ];

% check dimensions and type of image and edge map
[h,w,~]=size(I); assert(isa(I,'uint8') && size(I,3)==3);
if(nargin<2 || isempty(E)), E=zeros(h,w,'single'); end
assert(isa(E,'single') && size(E,1)==h && size(E,2)==w);
I=rgbConvert(I,'rgb');

if( type=='w' )
  % run watershed algorithm
  S = uint32(watershed(convTri(E,1))); b=1;
  
else
  if( ~isempty(o.seed) )
    % utilize seed segmentation removing boundaries if necessary
    S = o.seed; assert(isa(S,'uint32') && size(S,1)==h && size(S,2)==w);
    if(o.bounds), S = spDetectMex('boundaries',S,E,0,o.nThreads); end
    
  else
    % initialize superpixels at half resolution
    s=1/2; h1 = h-mod(h,1/s); w1 = w-mod(w,1/s);
    I0 = imResample(I(1:h1,1:w1,:),s);
    E0 = imResample(E(1:h1,1:w1),s);
    S = uint32(reshape(0:h1*w1*s*s-1,h1*s,w1*s));
    
    % refine superpixels at half resolution
    p = [o.nIter*2 o.nThreads sigs(1)*s*s sigs(2)/s/s sigs(3:4)];
    S = spDetectMex('sticky',S,convTri(I0,1),E0,p);
    S = imResample(S,1/s,'nearest');
    S = uint32(imPad(single(S),[0 h-h1 0 w-w1],'replicate'));
  end
  
  % refine superpixels at full resolution
  p = [o.nIter o.nThreads sigs]; b=0;
  S = spDetectMex('sticky',S,convTri(I,1),E,p);
  
end

% add or remove superpixel boundaries as necessary
if(o.bounds~=b), S = spDetectMex('boundaries',S,E,o.bounds,o.nThreads); end

% optionally merge superpixels
if(o.merge>0 && o.bounds), S = spDetectMex('merge',S,E,o.merge); end

% optionally create visualization
if(nargout>=2), V=spDetectMex('visualize',S,I,o.bounds); end

end
