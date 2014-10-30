function varargout = edgesEvalDir( varargin )
% Calculate edge precision/recall results for directory of edge images.
%
% Enhanced replacement for boundaryBench() from BSDS500 code:
%  http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/
% Uses same format for results and is fully compatible with boundaryBench.
% Given default prms results are *identical* to boundaryBench with the
% additional 9th output of R50 (recall at 50% precision).
%
% The odsF/P/R/T are results at the ODS (optimal dataset scale).
% The oisF/P/R are results at the OIS (optimal image scale).
% Naming convention: P=precision, R=recall, F=2/(1/P+1/R), T=threshold.
%
% In addition to the outputs, edgesEvalDir() creates three files:
%  eval_bdry_img.txt - per image OIS results [imgId T R P F]
%  eval_bdry_thr.txt - per threshold ODS results [T R P F]
%  eval_bdry.txt     - complete results (*re-ordered* copy of output)
% These files are identical to the ones created by boundaryBench.
%
% USAGE
%  [odsF,odsP,odsR,odsT,oisF,oisP,oisR,AP,R50] = edgesEvalDir( prms )
%  [ODS,~,~,~,OIS,~,~,AP,R50] = edgesEvalDir( prms )
%
% INPUTS
%  prms       - parameters (struct or name/value pairs)
%   .resDir     - ['REQ'] dir containing edge detection results (.png)
%   .gtDir      - ['REQ'] dir containing ground truth (.mat)
%   .pDistr     - [{'type','parfor'}] parameters for fevalDistr
%   .cleanup    - [0] if true delete temporary files
%   .thrs       - [99] number or vector of thresholds for evaluation
%   .maxDist    - [.0075] maximum tolerance for edge match
%   .thin       - [1] if true thin boundary maps
%
% OUTPUTS
%   odsF/P/R/T  - F-measure, precision, recall and threshold at ODS
%   oisF/P/R    - F-measure, precision, and recall at OIS
%   AP          - average precision
%   R50         - recall at 50% precision
%
% EXAMPLE
%
% See also edgesEvalImg, edgesEvalPlot
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

% get additional parameters
dfs={ 'resDir','REQ', 'gtDir','REQ', 'pDistr',{{'type','parfor'}}, ...
  'cleanup',0, 'thrs',99, 'maxDist',.0075, 'thin',1 };
p=getPrmDflt(varargin,dfs,1); resDir=p.resDir; gtDir=p.gtDir;
evalDir=[resDir '-eval/']; if(~exist(evalDir,'dir')), mkdir(evalDir); end

% check if results already exist, if so load and return
fNm = fullfile(evalDir,'eval_bdry.txt');
if(exist(fNm,'file')), R=dlmread(fNm); R=mat2cell2(R,[1 8]);
  varargout=R([4 3 2 1 7 6 5 8]); if(nargout<=8), return; end;
  R=dlmread(fullfile(evalDir,'eval_bdry_thr.txt')); P=R(:,3); R=R(:,2);
  [~,o]=unique(P); R50=interp1(P(o),R(o),max(P(o(1)),.5));
  varargout=[varargout R50]; return;
end

% perform evaluation on each image (this part can be very slow)
assert(exist(resDir,'dir')==7); assert(exist(gtDir,'dir')==7);
ids=dir(fullfile(gtDir,'*.mat')); ids={ids.name}; n=length(ids);
do=false(1,n); jobs=cell(1,n); res=cell(1,n);
for i=1:n, id=ids{i}(1:end-4);
  res{i}=fullfile(evalDir,[id '_ev1.txt']); do(i)=~exist(res{i},'file');
  im1=fullfile(resDir,[id '.png']); gt1=fullfile(gtDir,[id '.mat']);
  jobs{i}={im1,gt1,'out',res{i},'thrs',p.thrs,'maxDist',p.maxDist,...
    'thin',p.thin}; if(0), edgesEvalImg(jobs{i}{:}); end
end
fevalDistr('edgesEvalImg',jobs(do),p.pDistr{:});

% collect evaluation results
I=dlmread(res{1}); T=I(:,1);
Z=zeros(numel(T),1); cntR=Z; sumR=Z; cntP=Z; sumP=Z;
oisCntR=0; oisSumR=0; oisCntP=0; oisSumP=0; scores=zeros(n,5);
for i=1:n
  % load image results and accumulate
  I = dlmread(res{i});
  cntR1=I(:,2); cntR=cntR+cntR1; sumR1=I(:,3); sumR=sumR+sumR1;
  cntP1=I(:,4); cntP=cntP+cntP1; sumP1=I(:,5); sumP=sumP+sumP1;
  % compute OIS scores for image
  [R,P,F] = computeRPF(cntR1,sumR1,cntP1,sumP1); [~,k]=max(F);
  [oisR1,oisP1,oisF1,oisT1] = findBestRPF(T,R,P);
  scores(i,:) = [i oisT1 oisR1 oisP1 oisF1];
  % oisCnt/Sum will be used to compute dataset OIS scores
  oisCntR=oisCntR+cntR1(k); oisSumR=oisSumR+sumR1(k);
  oisCntP=oisCntP+cntP1(k); oisSumP=oisSumP+sumP1(k);
end

% compute ODS R/P/F and OIS R/P/F
[R,P,F] = computeRPF(cntR,sumR,cntP,sumP);
[odsR,odsP,odsF,odsT] = findBestRPF(T,R,P);
[oisR,oisP,oisF] = computeRPF(oisCntR,oisSumR,oisCntP,oisSumP);

% compute AP/R50 (interpolating 100 values, has minor bug: should be /101)
if(0), R=[0; R; 1]; P=[1; P; 0]; F=[0; F; 0]; T=[1; T; 0]; end
[~,k]=unique(R); k=k(end:-1:1); R=R(k); P=P(k); T=T(k); F=F(k); AP=0;
if(numel(R)>1), AP=interp1(R,P,0:.01:1); AP=sum(AP(~isnan(AP)))/100; end
[~,o]=unique(P); R50=interp1(P(o),R(o),max(P(o(1)),.5));

% write results to disk
varargout = {odsF,odsP,odsR,odsT,oisF,oisP,oisR,AP,R50};
writeRes(evalDir,'eval_bdry_img.txt',scores);
writeRes(evalDir,'eval_bdry_thr.txt',[T R P F]);
writeRes(evalDir,'eval_bdry.txt',[varargout{[4 3 2 1 7 6 5 8]}]);

% optionally perform cleanup
if( p.cleanup ), delete([evalDir '/*_ev1.txt']);
  delete([resDir '/*.png']); rmdir(resDir); end

end

function [R,P,F] = computeRPF(cntR,sumR,cntP,sumP)
% compute precision, recall and F measure given cnts and sums
R=cntR./max(eps,sumR); P=cntP./max(eps,sumP); F=2*P.*R./max(eps,P+R);
end

function [bstR,bstP,bstF,bstT] = findBestRPF(T,R,P)
% linearly interpolate to find best thr for optimizing F
if(numel(T)==1), bstT=T; bstR=R; bstP=P;
  bstF=2*P.*R./max(eps,P+R); return; end
A=linspace(0,1,100); B=1-A; bstF=-1;
for j = 2:numel(T)
  Rj=R(j).*A+R(j-1).*B; Pj=P(j).*A+P(j-1).*B; Tj=T(j).*A+T(j-1).*B;
  Fj=2.*Pj.*Rj./max(eps,Pj+Rj); [f,k]=max(Fj);
  if(f>bstF), bstT=Tj(k); bstR=Rj(k); bstP=Pj(k); bstF=f; end
end
end

function writeRes( alg, fNm, vals )
% write results to disk
k=size(vals,2); fNm=fullfile(alg,fNm); fid=fopen(fNm,'w');
if(fid==-1), error('Could not open file %s for writing.',fNm); end
frmt=repmat('%10g ',[1 k]); frmt=[frmt(1:end-1) '\n'];
fprintf(fid,frmt,vals'); fclose(fid);
end
