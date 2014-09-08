function varargout = edgesEval( model, varargin )
% Run and evaluate structured edge detector on BSDS500.
%
% This function first runs the trained structured edge detector on every
% test or validation image in BSDS then call edgesEvalDir.m to perform the
% actual edge evaluation. edgesEval is specific to the structured edge
% detector and BSDS, edgesEvalDir is general purpose edge evaluation code.
% For example usage of edgesEval see edgesDemo.
%
% USAGE
%  varargout = edgesEval( model, prms )
%
% INPUTS
%  model      - structured edge model trained with edgesTrain
%  prms       - parameters (struct or name/value pairs)
%   .dataType   - ['test'] should be either 'test' or 'val'
%   .name       - [''] name to append to evaluation
%   .opts       - {} list of model opts to overwrite
%   .show       - [0] if true plot results using edgesEvalPlot
%   .pDistr     - [{'type','parfor'}] parameters for fevalDistr
%   .cleanup    - [0] if true delete temporary files
%   .thrs       - [99] number or vector of thresholds for evaluation
%   .maxDist    - [.0075] maximum tolerance for edge match
%
% OUTPUTS
%  varargout  - same outputs as edgesEvalDir
%
% EXAMPLE
%
% See also edgesDemo, edgesDetect, edgesEvalDir, edgesEvalPlot
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

% get default parameters
dfs={'dataType','test', 'name','', 'opts',{}, 'show',0, ...
  'pDistr',{{'type','parfor'}}, 'cleanup',0, 'thrs',99, 'maxDist',.0075 };
p=getPrmDflt(varargin,dfs,1);

% load model and update model.opts acoording to opts
if( ischar(model) ), model=load(model); model=model.model; end
for i=1:length(p.opts)/2, model.opts.(p.opts{i*2-1})=p.opts{i*2}; end
rgbd=model.opts.rgbd; model.opts.nms=1;

% get list of relevant directories (image, depth, gt, results)
name = [model.opts.modelFnm,p.name];
imgDir = fullfile(model.opts.bsdsDir,'images',p.dataType);
depDir = fullfile(model.opts.bsdsDir,'depth',p.dataType);
gtDir  = fullfile(model.opts.bsdsDir,'groundTruth',p.dataType);
resDir = fullfile(model.opts.modelDir,p.dataType,name);
assert(exist(imgDir,'dir')==7); assert(exist(gtDir,'dir')==7);

% run edgesDetect() on every image in imgDir and store in resDir
if( ~exist(fullfile([resDir '-eval'],'eval_bdry.txt'),'file') )
  if(~exist(resDir,'dir')), mkdir(resDir); end
  ids=dir(imgDir); ids=ids([ids.bytes]>0); ids={ids.name}; n=length(ids);
  ext=ids{1}(end-2:end); for i=1:n, ids{i}=ids{i}(1:end-4); end
  res=cell(1,n); for i=1:n, res{i}=fullfile(resDir,[ids{i} '.png']); end
  do=false(1,n); for i=1:n, do(i)=~exist(res{i},'file'); end
  ids=ids(do); res=res(do); m=length(ids);
  parfor i=1:m, id=ids{i};
    I = imread(fullfile(imgDir,[id '.' ext])); D=[];
    if(rgbd), D=single(imread(fullfile(depDir,[id '.png'])))/1e4; end
    if(rgbd==1), I=D; elseif(rgbd==2), I=cat(3,single(I)/255,D); end
    E=edgesDetect(I,model); imwrite(uint8(E*255),res{i});
  end
end

% perform actual evaluation using edgesEvalDir
varargout=cell(1,max(1,nargout));
[varargout{:}] = edgesEvalDir('resDir',resDir,'gtDir',gtDir,...
  'pDistr',p.pDistr,'cleanup',p.cleanup,'thrs',p.thrs,'maxDist',p.maxDist);
if( p.show ), figure(p.show); edgesEvalPlot(resDir,name); end

end
