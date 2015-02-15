function edgeBoxesSweeps()
% Parameter sweeps for Edges Boxes object proposals.
%
% Running the parameter sweeps requires altering internal flags.
% The sweeps are not well documented, use at your own discretion.
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar and Larry Zitnick, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

% define parameter sweeps
rt = [fileparts(mfilename('fullpath')) filesep];
expNms = {'alpha','beta','eta','minScore','edgeMinMag','edgeMergeThr',...
  'clusterMinMag','maxAspectRatio','minBoxArea','gamma','kappa'};
expNms=expNms(1:end); opts=createExp(rt,expNms); maxn=inf;

% run training and testing jobs
jobs = createJobs(rt,opts,maxn);
tic, for i=1:length(jobs), edgeBoxes(jobs{i}{:}); end; toc

% plot all results
for e=1:length(expNms)
  eval={'data',boxesData('split','val'),'names',{opts{e}.name},...
    'resDir',[rt 'boxes/sweeps/'],'maxn',maxn,'fName',expNms{e}};
  boxesEval(eval{:},'thrs',.7);
  ar=boxesEval(eval{:},'thrs',.5:.05:1,'cnts',1000);
  ar=squeeze(mean(ar,2)); [~,i]=max(ar); disp(opts{e}(i))
end

end

function jobs = createJobs( rt, opts, maxn )
% create jobs
M='models/forest/modelBsds'; M=load(M); M=M.model; M.opts.nThreads=1;
data=boxesData('split','val'); fs=data.imgs(1:min(end,maxn));
opts=[opts{:}]; N=length(opts); jobs=cell(1,N); D=zeros(1,N);
for e=1:N, opts(e).name=[rt 'boxes/sweeps/' opts(e).name '-val.mat']; end
for e=1:N, D(e)=exist(opts(e).name,'file')==2; jobs{e}={fs,M,opts(e)}; end
[~,K]=unique({opts.name},'stable'); D=D(K); jobs=jobs(K); jobs=jobs(~D);
fprintf('nJobs = %i\n',length(jobs));
end

function opts = createExp( rt, expNm )

% if expNm is a cell, call recursively and return
if( iscell(expNm) )
  N=length(expNm); opts=cell(1,N);
  for e=1:N, opts{e}=createExp(rt,expNm{e}); end; return;
end

% setup opts
opts=edgeBoxes(); opts.minScore=0;
N=100; optsDefault=opts; opts=opts(ones(1,N));
switch expNm
  case 'alpha'
    vs=45:5:75; N=length(vs);
    for e=1:N, opts(e).alpha=vs(e)/100; end
  case 'beta'
    vs=60:5:90; N=length(vs);
    for e=1:N, opts(e).beta=vs(e)/100; end
  case 'eta'
    vs=0:5; N=length(vs);
    for e=1:N, opts(e).eta=1-vs(e)/10000; end
  case 'minScore'
    vs=[0 5 10 25 50 100]; N=length(vs);
    for e=1:N, opts(e).minScore=vs(e)/1000; end
  case 'edgeMinMag'
    vs=[0 50 100 200 400]; N=length(vs);
    for e=1:N, opts(e).edgeMinMag=vs(e)/1000; end
  case 'edgeMergeThr'
    vs=[25 50 100 200 400]; N=length(vs);
    for e=1:N, opts(e).edgeMergeThr=vs(e)/100; end
  case 'clusterMinMag'
    vs=[0 50 100 200 400]; N=length(vs);
    for e=1:N, opts(e).clusterMinMag=vs(e)/100; end
  case 'maxAspectRatio'
    vs=1:5; N=length(vs);
    for e=1:N, opts(e).maxAspectRatio=vs(e); end
  case 'minBoxArea'
    vs=[100 250 500 1000 2500 5000]; N=length(vs);
    for e=1:N, opts(e).minBoxArea=vs(e); end
  case 'gamma'
    vs=[25 50 100 200 400 10000]; N=length(vs);
    for e=1:N, opts(e).gamma=vs(e)/100; end
  case 'kappa'
    vs=50:25:200; N=length(vs);
    for e=1:N, opts(e).kappa=vs(e)/100; end
  otherwise, error('invalid exp: %s',expNm);
end

% produce final set of opts and find default opts
O=1:N; opts=opts(O); d=0;
for e=1:N, if(isequal(optsDefault,opts(e))), d=e; break; end; end
if(d==0), disp(expNm); assert(false); end; opts(d).name='Default';
for e=1:N, if(e~=d), opts(e).name=[expNm int2str2(vs(e),5)]; end; end

end
