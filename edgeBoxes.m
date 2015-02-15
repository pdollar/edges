function bbs = edgeBoxes( I, model, varargin )
% Generate Edge Boxes object proposals in given image(s).
%
% Compute Edge Boxes object proposals as described in:
%  C. Lawrence Zitnick and Piotr Dollár
%  "Edge Boxes: Locating Object Proposals from Edges", ECCV 2014.
% The proposal boxes are fast to compute and give state-of-the-art recall.
% Please cite the above paper if you end up using the code.
%
% The most important params are alpha and beta. The defaults are optimized
% for detecting boxes at intersection over union (IoU) of 0.7. For other
% settings of alpha/beta see the ECCV paper. In general larger alpha/beta
% improve results at higher IoU (but using large alpha can be quite slow).
% minScore/maxBoxes control the number of boxes returned and impact speed.
% Finally, a number of additional params listed below are set to reasonable
% defaults and in most cases should not need to be altered.
%
% We recently introduced an adaptive variant of nms that results in better
% Average Recall (AR) and also better subsequent detection performance.
% This variant is described in Section 5.E of the following paper:
%  Jan Hosang, Rodrigo Benenson, Piotr Dollár, and Bernt Schiele
%  "What makes for effective detection proposals?", arXiv 2015.
% TL;DR: to get top AR for 1000 boxes set alpha=.65, beta=.90, eta=.9996.
%
% For a faster version the proposal code runs at ~10 fps on average use:
%  model.opts.sharpen=0; opts.alpha=.625; opts.minScore=.02;
%
% The code uses the Structured Edge Detector to compute edge strength and
% orientation. See edgesDetect.m for details. Alternatively, the code could
% be altered to use any other edge detector such as Canny.
%
% The input 'I' can either be a single (color) image (or filename) or a
% cell array of images (or filenames). In the first case, the return is a
% set of bbs where each row has the format [x y w h score] and score is the
% confidence of detection. If the input is a cell array, the output is a
% cell array where each element is a set of bbs in the form above (in this
% case a parfor loop is used to speed execution).
%
% USAGE
%  opts = edgeBoxes()
%  bbs = edgeBoxes( I, model, opts )
%
% INPUTS
%  I          - input image(s) of filename(s) of input image(s)
%  model      - Structured Edge model trained with edgesTrain
%  opts       - parameters (struct or name/value pairs)
%   (1) main parameters, see above for details
%   .name           - [] target filename (if specified return is 1)
%   .alpha          - [.65] step size of sliding window search
%   .beta           - [.75] nms threshold for object proposals
%   .eta            - [1.0] adaptation rate for nms threshold (see arXiv15)
%   .minScore       - [.01] min score of boxes to detect
%   .maxBoxes       - [1e4] max number of boxes to detect
%   (2) additional parameters, safe to ignore and leave at default vals
%   .edgeMinMag     - [.1] increase to trade off accuracy for speed
%   .edgeMergeThr   - [.5] increase to trade off accuracy for speed
%   .clusterMinMag  - [.5] increase to trade off accuracy for speed
%   .maxAspectRatio - [3] max aspect ratio of boxes
%   .minBoxArea     - [1000] minimum area of boxes
%   .gamma          - [2] affinity sensitivity, see equation 1 in paper
%   .kappa          - [1.5] scale sensitivity, see equation 3 in paper
%
% OUTPUTS
%  bbs        - [nx5] array containing proposal bbs [x y w h score]
%
% EXAMPLE
%
% See also edgeBoxesDemo, edgesDetect
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar and Larry Zitnick, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

% get default parameters (unimportant parameters are undocumented)
dfs={'name','', 'alpha',.65, 'beta',.75, 'eta',1, 'minScore',.01, ...
  'maxBoxes',1e4, 'edgeMinMag',.1, 'edgeMergeThr',.5,'clusterMinMag',.5,...
  'maxAspectRatio',3, 'minBoxArea',1000, 'gamma',2, 'kappa',1.5 };
o=getPrmDflt(varargin,dfs,1); if(nargin==0), bbs=o; return; end

% run detector possibly over multiple images and optionally save results
f=o.name; if(~isempty(f) && exist(f,'file')), bbs=1; return; end
if(~iscell(I)), bbs=edgeBoxesImg(I,model,o); else n=length(I);
  bbs=cell(n,1); parfor i=1:n, bbs{i}=edgeBoxesImg(I{i},model,o); end; end
d=fileparts(f); if(~isempty(d)&&~exist(d,'dir')), mkdir(d); end
if(~isempty(f)), save(f,'bbs'); bbs=1; end

end

function bbs = edgeBoxesImg( I, model, o )
% Generate Edge Boxes object proposals in single image.
if(all(ischar(I))), I=imread(I); end
model.opts.nms=0; [E,O]=edgesDetect(I,model);
if(0), E=gradientMag(convTri(single(I),4)); E=E/max(E(:)); end
E=edgesNmsMex(E,O,2,0,1,model.opts.nThreads);
bbs=edgeBoxesMex(E,O,o.alpha,o.beta,o.eta,o.minScore,o.maxBoxes,...
  o.edgeMinMag,o.edgeMergeThr,o.clusterMinMag,...
  o.maxAspectRatio,o.minBoxArea,o.gamma,o.kappa);
end
