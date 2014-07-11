function recall = boxesEval( varargin )
% Perform object proposal bounding box evaluation.

cs=cell(100,1); for i=1:100, cs{i}=max(.3,mod([.3 .47 .16]*(i+1),1)); end
ws=[1 2 5 10 20 50 100 200 500 1000 2000 5000];
dfs={ 'detectors',{}, 'split','val', 'maxn',inf, 'fName','', 'show',1, ...
  'thrs',.5:.2:.9, 'windows',ws, 'clrs',cs, 'resDir','boxes/' };
o=getPrmDflt(varargin,dfs,1); if(nargin==0), recall=o; return; end
recall=boxesEvalAll(o); if(o.show), plotResult(recall,o); end

end

function recall = boxesEvalAll( o )
% compute and gather all results (caches individual results to disk)
W=length(o.windows); T=length(o.thrs); K=length(o.detectors);
gt=boxesData('resDir',o.resDir,'split',o.split); n=min(o.maxn,gt.n);
gt=gt.gt(1:n); recall=zeros(W,T,K); [ws,ts,ks]=ndgrid(1:W,1:T,1:K);
parfor i=1:W*T*K, w=ws(i); t=ts(i); k=ks(i);
  % if evaluation result exists simply load it
  rdir=[o.resDir '/eval/' o.detectors{k} '/' o.split '/'];
  rnm=[rdir 'N' int2str2(n,5) '-W' int2str2(o.windows(w),5) ...
    '-T' int2str2(round(o.thrs(t)*100),2) '.txt']; %#ok<*PFBNS>
  if(exist(rnm,'file')), recall(i)=load(rnm,'-ascii'); continue; end
  % perform evaluation if result does not exist
  dt=load([o.resDir '/' o.detectors{k} '-' o.split]);
  if(isfield(dt,'dt')), dt=dt.dt(1:n); else dt=dt.bbs(1:n); end
  dt1=dt; for j=1:n, dt1{j}=dt1{j}(1:min(end,o.windows(w)),:); end
  [gt1,dt1]=bbGt('evalRes',gt,dt1,o.thrs(t));
  [~,r]=bbGt('compRoc',gt1,dt1,1); r=max(r); recall(i)=r;
  if(~exist(rdir,'dir')), mkdir(rdir); end; dlmwrite(rnm,r);
end
% display summary statistics
[ts,ks]=ndgrid(1:T,1:K); ms=log(o.windows); rt=.75;
for i=1:T*K, t=ts(i); k=ks(i); r=recall(:,t,k)'; if(W==1), continue; end
  a=find(rt<=r); if(isempty(a)), m=inf; else a=a(1); b=a-1;
    m=round(exp((rt-r(b))/(r(a)-r(b))*(ms(a)-ms(b))+ms(b))); end
  auc=sum(diff(ms/ms(end)).*(r(1:end-1)+r(2:end))/2);
  fprintf('%15s  T=%.2f  A=%.2f  W=%4i  R=%.2f\n',...
    o.detectors{k},o.thrs(t),auc,m,max(r));
end
% optionally save results to text file
if(isempty(o.fName)), return; end
d=[o.resDir '/plots/']; if(~exist(d,'dir')), mkdir(d); end
dlmwrite([d o.fName '-' o.split '.txt'],squeeze(recall));
end

function plotResult( recall, o )
% plot results
[W,T,K]=size(recall); fSiz={'FontSize',12};
for type=1:2
  if(type==1), xs=o.windows; else xs=o.thrs; end;
  if(length(xs)==1), continue; end; s=[T,W]; M=s(type);
  R=recall; if(type==2), R=permute(R,[2 1 3]); end
  figure(o.show+type-1); clf; hold on; hs=zeros(M,K);
  for i=1:M, for k=1:K, hs(i,k)=plot(xs,R(:,i,k),...
        'Color',o.clrs{k},'LineWidth',3); end; end
  s={'# of windows','IoU'}; xlabel(s{type},fSiz{:});
  s={'log','linear'}; set(gca,'XScale',s{type});
  ylabel('Detection Rate',fSiz{:}); set(gca,'YTick',0:.2:1);
  hold off; axis([min(xs) max(xs) 0 1]); grid on; set(gca,fSiz{:});
  set(gca,'XMinorGrid','off','XMinorTic','off');
  set(gca,'YMinorGrid','off','YMinorTic','off');
  s={'nw','ne'}; legend(hs(1,:),o.detectors,'Location',s{type});
  if(isempty(o.fName)), continue; end; s={'Win','IoU'};
  savefig([o.resDir '/plots/' s{type} '-' o.split '-' o.fName],'png');
end
end
