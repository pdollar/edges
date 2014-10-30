function edgesEvalPlot( algs, nms, cols )
% Plot edge precision/recall results for directory of edge images.
%
% Enhanced replacement for plot_eval() from BSDS500 code:
%  http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/
% Uses same format and is fully compatible with plot_eval. Use this
% function to plot the edge results created using edgesEvalDir.
%
% USAGE
%  edgesEvalPlot( algs, [nms], [cols] )
%
% INPUTS
%  algs       - {nx1} algorithm result directories
%  nms        - [{nx1}] algorithm names (for legend)
%  cols       - [{nx1}] algorithm colors
%
% OUTPUTS
%
% EXAMPLE
%
% See also edgesEvalDir
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

% parse inputs
if(nargin<2||isempty(nms)), nms={}; end; if(~iscell(nms)), nms={nms}; end
if(nargin<3||isempty(cols)), cols=repmat({'r','g','b','k','m'},1,100); end
if(~iscell(algs)), algs={algs}; end; if(~iscell(cols)), cols={cols}; end

% setup basic plot (isometric contour lines and human performance)
clf; box on; grid on; hold on;
line([0 1],[.5 .5],'LineWidth',2,'Color',.7*[1 1 1]);
for f=0.1:0.1:0.9, r=f:0.01:1; p=f.*r./(2.*r-f); %f=2./(1./p+1./r)
  plot(r,p,'Color',[0 1 0]); plot(p,r,'Color',[0 1 0]); end
if(1), h=plot(0.7235,0.9014,'o','MarkerSize',8,'Color',[0 .5 0],...
    'MarkerFaceColor',[0 .5 0],'MarkerEdgeColor',[0 .5 0]); end
set(gca,'XTick',0:0.1:1,'YTick',0:0.1:1);
grid on; xlabel('Recall'); ylabel('Precision');
axis equal; axis([0 1 0 1]);

% load results for every algorithm (pr=[T,R,P,F])
n=length(algs); hs=zeros(1,n); res=zeros(n,9); prs=cell(1,n);
for i=1:n, a=[algs{i} '-eval'];
  pr=dlmread(fullfile(a,'eval_bdry_thr.txt')); pr=pr(pr(:,2)>=1e-3,:);
  [~,o]=unique(pr(:,3)); R50=interp1(pr(o,3),pr(o,2),max(pr(o(1),3),.5));
  res(i,1:8)=dlmread(fullfile(a,'eval_bdry.txt')); res(i,9)=R50; prs{i}=pr;
end;

% sort algorithms by ODS score
[~,o]=sort(res(:,4),'descend'); res=res(o,:); prs=prs(o);
cols=cols(o); if(~isempty(nms)), nms=nms(o); end

% plot results for every algorithm (plot best last)
for i=n:-1:1
  hs(i)=plot(prs{i}(:,2),prs{i}(:,3),'-','LineWidth',3,'Color',cols{i});
  fprintf('ODS=%.3f OIS=%.3f AP=%.3f R50=%.3f',res(i,[4 7:9]));
  if(~isempty(nms)), fprintf(' - %s',nms{i}); end; fprintf('\n');
end

% show legend if nms provided (report best first)
hold off; if(isempty(nms)), return; end
for i=1:n, nms{i}=sprintf('[F=.%i] %s',round(res(i,4)*100),nms{i}); end
if(1), hs=[h hs]; nms=['[F=.80] Human'; nms(:)]; end
legend(hs,nms,'Location','sw');

end
