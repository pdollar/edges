function edgesSweeps()
% Parameter sweeps for structured edge detector.
%
% Running the parameter sweeps requires altering internal flags.
% The sweeps are not well documented, use at your own discretion.
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

% select type and location of cluster (see fevalDistr.m)
rtDir = 'D:\code\research\edges\';
pDistrTrn={'type','local'}; pDistrTst=pDistrTrn;

% define parameter sweeps
expNms= {'MD-imWidth','MD-gtWidth','TR-nData','TR-nPos','TR-nImgs',...
  'TR-sharpen','TR-nTrees','TR-split','TR-minChild','TR-maxDepth',...
  'TR-fracFtrs','TR-nSamples','TR-nClasses','TR-discretize',...
  'FT-nCells','FT-nOrients','FT-normRad','FT-shrink',...
  'FT-grdSmooth','FT-chnSmooth','FT-simSmooth','final'};
expNms=expNms(1:end); T=5; full=3; dataset='BSDS';
[opts,lgd,lbl]=createExp(rtDir,expNms,dataset,full);

% run training and testing jobs
[jobsTrn,jobsTst] = createJobs(rtDir,opts,dataset,T); N=length(expNms);
fprintf('nTrain = %i; nTest = %i\n',length(jobsTrn),length(jobsTst));
tic, s=fevalDistr('edgesTrain',jobsTrn,pDistrTrn); assert(s==1); toc
tic, s=fevalDistr('edgesEval',jobsTst,pDistrTst); assert(s==1); toc

% plot results
for e=1:N, plotExps(expNms{e},opts{e},lgd{e},lbl{e},T,1); end

end

function plotExps( expNm, opts, lgd, lbl, T, type )
% get all results and display error
disp([expNm ' [' lbl ']']); N=length(lgd);
res=zeros(N,T); mNms=cell(1,N);
for e=1:N, mNms{e}=[opts(e).modelDir 'val/' opts(e).modelFnm]; end
for e=1:N, for t=1:T, r=dlmread([mNms{e} 'T' int2str2(t,2) ...
      '-eval/eval_bdry.txt']); r=r([4 7 8]); res(e,t)=r(type); end; end
stds=std(res,0,2)*100; R=mean(res,2)*100; msg=' %.2f +/- %.2f  [%s]\n';
for e=1:N, fprintf(msg,R(e),stds(e),lgd{e}); end
if(0), disp(res); disp(max(res,[],2)); end
types={'ODS','OIS','AP'}; type=types{type};
% plot sweeps (two cases for format of x labels)
figPrp = {'Units','Pixels','Position',[800 600 640 220]};
figure(1); clf; set(1,figPrp{:}); set(gca,'FontSize',24); clr=[0 .69 .94];
pPl1={'LineWidth',3,'MarkerSize',15,'Color',clr,'MarkerFaceColor',clr};
pPl2=pPl1; clr=[1 .75 0]; pPl2{6}=clr; pPl2{8}=clr; d=0;
for e=1:N, if(lgd{e}(end)=='*'), d=e; end; end; if(d), lgd{d}(end)=[]; end
plot(R,'-d',pPl1{:}); hold on; if(d),plot(d,R(d),'d',pPl2{:}); end; e=.001;
ylabel([type ' \times 100']); axis([.5 N+.5 min([R; 66])-e max([R; 72])+e])
if(isempty(lbl)), imLabel(lgd,'bottom',30,{'FontSize',24}); lgd=[]; end
if(0); xlabel(lbl); end; set(gca,'XTick',1:1:N,'XTickLabel',lgd(1:1:N));
% save plot
plDir=[opts(1).modelDir 'plots' type '/']; fFig=[plDir expNm];
if(~exist(plDir,'dir')), mkdir(plDir); end %#ok<*CTCH>
for t=1:25, try savefig(fFig,1,'png'); break; catch, pause(1), end; end
end

function [jobsTrn,jobsTst] = createJobs( rtDir, opts, dataset, T )
% Prepare all jobs (one train and one test job per set of opts).
opts=[opts{:}]; N=length(opts); NT=N*T;
opts=repmat(opts,1,T); nms=cell(1,NT);
jobsTrn=cell(1,NT); doneTrn=zeros(1,NT);
jobsTst=cell(1,NT); doneTst=zeros(1,NT);
thrs = logspace(-2,log10(.99),25)'; thrs(1)=1e-5;
if( strcmpi(dataset,'bsds') )
  pTest={'dataType','val', 'thrs',thrs, 'cleanup',1,...
    'opts',{'modelDir',[rtDir '/sweepsBSDS/'],...
    'bsdsDir',[rtDir '/BSR/BSDS500/data/']} };
else
  pTest={'dataType','val', 'thrs',thrs, 'cleanup',1, 'maxDist',.011,...
    'opts',{'modelDir',[rtDir '/sweepsNYUD/'],...
    'bsdsDir',[rtDir '/BSR/NYUD/data/']} };
end
for e=1:NT
  t=ceil(e/N); opts(e).seed=(t-1)*100000+1;
  nm=[opts(e).modelFnm 'T' int2str2(t,2)]; opts(e).modelFnm=nm;
  mFnm=[opts(e).modelDir 'forest/' nm '.mat']; nms{e}=nm;
  eFnm=[opts(e).modelDir 'val/' nm '-eval/eval_bdry.txt'];
  doneTrn(e)=exist(mFnm,'file')==2; jobsTrn{e}={opts(e)};
  doneTst(e)=exist(eFnm,'file')==2; jobsTst{e}=[mFnm pTest];
end
[~,kp]=unique(nms,'stable');
doneTrn=doneTrn(kp); jobsTrn=jobsTrn(kp); jobsTrn=jobsTrn(~doneTrn);
doneTst=doneTst(kp); jobsTst=jobsTst(kp); jobsTst=jobsTst(~doneTst);
end

function [opts,lgd,lbl] = createExp( rtDir, expNm, dataset, full )

% if expNm is a cell, call recursively and return
if( iscell(expNm) )
  N=length(expNm); opts=cell(1,N); lgd=opts; lbl=opts;
  for e=1:N, [opts{e},lgd{e},lbl{e}]=...
      createExp(rtDir,expNm{e},dataset,full); end; return;
end

% default params for edgesTrain.m
opts=edgesTrain(); opts.nThreads=1; nData=2e5; nPos=.5;
opts.nPos=round(nData*nPos); opts.nNeg=round(nData*(1-nPos));
assert(any(strcmpi(dataset,{'bsds','nyud'})));
if( strcmpi(dataset,'bsds') )
  opts.modelDir=[rtDir '/sweepsBSDS/']; opts.nImgs=200;
  opts.bsdsDir=[rtDir '/BSR/BSDS500/data/'];
else
  opts.modelDir=[rtDir '/sweepsNYUD/']; opts.nImgs=381;
  opts.bsdsDir=[rtDir '/BSR/NYUD/data/'];
  opts.rgbd=2; opts.fracFtrs=1/8;
end

% setup opts
optsDefault=opts; N=100; lgd=cell(1,N); ss=lgd;
opts=opts(ones(1,N)); hasDefault=1;
switch expNm
  case 'MD-imWidth'
    lbl='window size for x'; vs=[1:4 6 8]*8; N=length(vs);
    for e=1:N, opts(e).imWidth=vs(e); end
    for e=1:N, opts(e).gtWidth=min(vs(e),opts(e).gtWidth); end
  case 'MD-gtWidth'
    lbl='window size for y'; vs=[1:4 6 8]*4; N=length(vs);
    for e=1:N, opts(e).gtWidth=vs(e); end
  case 'TR-nData'
    lbl='# train patches x10^4'; vs=[1 2 5 10 20 50 100 200]; N=length(vs);
    for e=1:N, opts(e).nPos=round(vs(e)*1e4*nPos); end
    for e=1:N, opts(e).nNeg=round(vs(e)*1e4*(1-nPos)); end
    if(full<2), N=5; elseif(full<3), N=6; end
  case 'TR-nPos'
    lbl = 'fraction positive data'; vs=20:10:80; N=length(vs);
    for e=1:N, opts(e).nPos=round(nData*vs(e)/100); end
    for e=1:N, opts(e).nNeg=round(nData*(1-vs(e)/100)); end
    for e=1:N, lgd{e}=sprintf('.%i',vs(e)/10); end
  case 'TR-nImgs'
    lbl='# train images'; vs=[10 20 50 100 opts(1).nImgs]; N=length(vs);
    for e=1:N, opts(e).nImgs=vs(e); end
  case 'TR-sharpen'
    lbl='sharpening radius'; vs=0:4; N=length(vs);
    for e=1:N, opts(e).sharpen=vs(e); end
  case 'TR-nTrees'
    lbl='# decision trees'; vs=2.^(0:4); N=length(vs);
    for e=1:N, opts(e).nTrees=vs(e); end;
    for e=1:N, opts(e).nTreesEval=max(1,vs(e)/2); end
    if(full<2), N=4; end
  case 'TR-split'
    lbl='information gain';
    ss={'gini','entropy','twoing'}; N=length(ss); lgd=ss;
    for e=1:N, opts(e).split=ss{e}; end
  case 'TR-minChild'
    lbl='min samples per node'; vs=2:2:16; N=length(vs);
    for e=1:N, opts(e).minChild=vs(e); end
  case 'TR-maxDepth'
    lbl='max tree depth'; vs=2.^(2:6); N=length(vs);
    for e=1:N, opts(e).maxDepth=vs(e); end
  case 'TR-fracFtrs'
    lbl='fraction features'; vs=2.^(2:6); N=length(vs);
    for e=1:N, opts(e).fracFtrs=1/vs(e); end
    for e=1:N, lgd{e}=sprintf('1/%i',vs(e)); end
  case 'TR-nSamples'
    lbl='m (size of Z)'; vs=2.^(0:2:8); N=length(vs);
    for e=1:N, opts(e).nSamples=vs(e); end
  case 'TR-nClasses'
    lbl='k (size of C)'; vs=2.^(1:5); N=length(vs);
    for e=1:N, opts(e).nClasses=vs(e); end
  case 'TR-discretize'
    lbl='discretization type';
    ss={'pca','kmeans'}; N=length(ss); lgd=ss;
    for e=1:N, opts(e).discretize=ss{e}; end
  case 'FT-nCells'
    lbl='# grid cells'; vs=1:2:7; N=length(vs);
    for e=1:N, opts(e).nCells=vs(e); end
    if(full<2), N=3; end
  case 'FT-nOrients'
    lbl='# gradient orients'; vs=0:2:8; N=length(vs);
    for e=1:N, opts(e).nOrients=vs(e); end
  case 'FT-normRad'
    lbl='normalization radius'; vs=[0 2.^(0:3)]; N=length(vs);
    for e=1:N, opts(e).normRad=vs(e); end
  case 'FT-shrink'
    lbl='channel downsample'; vs=2.^(0:2); N=length(vs);
    for e=1:N, opts(e).shrink=vs(e); end
  case 'FT-grdSmooth'
    lbl = 'gradient blur'; vs=[0 2.^(0:4)]; N=length(vs);
    for e=1:N, opts(e).grdSmooth=vs(e); end
  case 'FT-chnSmooth'
    lbl='channel blur'; vs=[0 2.^(0:4)]; N=length(vs);
    for e=1:N, opts(e).chnSmooth=vs(e); end
  case 'FT-simSmooth'
    lbl='self-similarity blur'; vs=[0 2.^(0:4)]; N=length(vs);
    for e=1:N, opts(e).simSmooth=vs(e); end
  case 'final'
    if(strcmpi(dataset,'bsds')), vs=[0 2]; else
      vs=0:3; rgbd=[2 2 1 0]; for e=2:4, opts(e).rgbd=rgbd(e); end; end
    lbl='final'; N=length(vs); nData=2e6;
    for e=2:N, opts(e).nPos=round(nData*nPos); end
    for e=2:N, opts(e).nNeg=round(nData*(1-nPos)); end
    if(full<3), N=1; end
  otherwise, error('invalid exp: %s',expNm);
end

% produce final set of opts and find default opts
for e=1:N, if(isempty(lgd{e})), lgd{e}=int2str(vs(e)); end; end
for e=1:N, if(isempty(ss{e})), ss{e}=int2str2(vs(e),5); end; end
O=1:N; opts=opts(O); lgd=lgd(O); ss=ss(O); d=0;
for e=1:N, if(isequal(optsDefault,opts(e))), d=e; break; end; end
if(hasDefault && d==0), disp(expNm); assert(false); end
for e=1:N, opts(e).modelFnm=[expNm ss{e}]; end
if(hasDefault), lgd{d}=[lgd{d} '*']; opts(d).modelFnm='Default'; end
if(0), disp([ss' lgd']'); end

end
