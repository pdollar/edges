% Demo for Sticky Superpixels (please see readme.txt first).

%% load pre-trained edge detection model and set opts (see edgesDemo.m)
model=load('models/forest/modelBsds'); model=model.model;
model.opts.nms=-1; model.opts.nThreads=4;
model.opts.multiscale=0; model.opts.sharpen=2;

%% set up opts for spDetect (see spDetect.m)
opts = spDetect;
opts.nThreads = 4;  % number of computation threads
opts.k = 512;       % controls scale of superpixels (big k -> big sp)
opts.alpha = .5;    % relative importance of regularity versus data terms
opts.beta = .9;     % relative importance of edge versus color terms
opts.merge = 0;     % set to small value to merge nearby superpixels at end

%% detect and display superpixels (see spDetect.m)
I = imread('peppers.png');
[E,~,~,segs]=edgesDetect(I,model);
tic, [S,V] = spDetect(I,E,opts); toc
figure(1); im(I); figure(2); im(V);

%% compute ultrametric contour map from superpixels (see spAffinities.m)
tic, [~,~,U]=spAffinities(S,E,segs,opts.nThreads); toc
figure(3); im(1-U); return;

%% compute video superpixels reusing initialization from previous frame
Is=seqIo(which('peds30.seq'),'toImgs'); Vs=single(Is); opts.bounds=0; tic
for i=1:size(Is,4), I=Is(:,:,:,i); E=edgesDetect(I,model);
  [opts.seed,Vs(:,:,:,i)]=spDetect(I,E,opts); end; opts.seed=[]; toc
Vs=uint8(Vs*255); playMovie([Is Vs],15,-10,struct('hasChn',1))
