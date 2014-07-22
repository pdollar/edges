% Demo for Sticky Superpixels (please see readme.txt first).

%% load pre-trained edge detection model and set opts (see edgesDemo.m)
model=load('models/forest/modelBsds'); model=model.model;
model.opts.multiscale=0; model.opts.sharpen=2; model.opts.nThreads=4;

%% set up opts for spDetect (see spDetect.m)
opts = spDetect;
opts.nThreads=4;  % number of computation threads
opts.k = 512;     % controls scale of superpixels (big k -> big sp)
opts.alpha = .5;  % relative importance of regularity versus data terms
opts.beta = .9;   % relative importance of edge versus color terms
opts.merge = 0;   % set to small value to merge nearby superpixels at end

%% detect and display superpixels (see spDetect.m)
I = imread('peppers.png');
[E,~,~,segs]=edgesDetect(I,model);
tic, [S,V] = spDetect(I,E,opts); toc
figure(1); im(V);

%% Convert superpixels boundaries back to edges (see spAffinities).
tic, [A,SE]=spAffinities(S,E,segs); toc
figure(2); im(1-SE);
