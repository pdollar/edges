% Demo for Edge Boxes (please see readme.txt first).

%% load pre-trained edge detection model (see edgesDemo.m)
model=load('models/forest/modelBsds'); model=model.model;

%% set up opts for edgeBoxes (see edgeBoxes.m)
opts = edgeBoxes;
opts.alpha = .65;     % step size of sliding window search
opts.beta  = .75;     % nms threshold for object proposals
opts.minScore = .01;  % min score of boxes to detect
opts.maxBoxes = 1e4;  % max number of boxes to detect

%% detect bbs (no visualization code for now)
I = imread('peppers.png');
tic, bbs=edgeBoxes(I,model,opts); toc
