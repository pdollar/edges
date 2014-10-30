% Demo for Edge Boxes (please see readme.txt first).

%% load pre-trained edge detection model and set opts (see edgesDemo.m)
model=load('models/forest/modelBsds'); model=model.model;
model.opts.multiscale=0; model.opts.sharpen=2; model.opts.nThreads=4;

%% set up opts for edgeBoxes (see edgeBoxes.m)
opts = edgeBoxes;
opts.alpha = .65;     % step size of sliding window search
opts.beta  = .75;     % nms threshold for object proposals
opts.minScore = .01;  % min score of boxes to detect
opts.maxBoxes = 1e4;  % max number of boxes to detect

%% detect Edge Box bounding box proposals (see edgeBoxes.m)
I = imread('peppers.png');
tic, bbs=edgeBoxes(I,model,opts); toc

%% show evaluation results (using pre-defined or interactive boxes)
gt=[122 248 92 65; 193 82 71 53; 410 237 101 81; 204 160 114 95; ...
  9 185 86 90; 389 93 120 117; 253 103 107 57; 81 140 91 63];
if(0), gt='Please select an object box.'; disp(gt); figure(1); imshow(I);
  title(gt); [~,gt]=imRectRot('rotate',0); gt=gt.getPos(); end
gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
figure(1); bbGt('showRes',I,gtRes,dtRes(dtRes(:,6)==1,:));
title('green=matched gt  red=missed gt  dashed-green=matched detect');

%% run and evaluate on entire dataset (see boxesData.m and boxesEval.m)
if(~exist('boxes/VOCdevkit/','dir')), return; end
split='val'; data=boxesData('split',split);
nm='EdgeBoxes70'; opts.name=['boxes/' nm '-' split '.mat'];
edgeBoxes(data.imgs,model,opts); opts.name=[];
boxesEval('data',data,'names',nm,'thrs',.7,'show',2);
boxesEval('data',data,'names',nm,'thrs',.5:.05:1,'cnts',1000,'show',3);
