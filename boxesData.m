function data = boxesData( varargin )
% Get ground truth data for object proposal bounding box evaluation.
%
% Used to get dataset information for evaluation of bounding box object
% proposals using boxesEval.m. This requires separate download of the
% appropriate dataset first. Currently only set up for PASCAL VOC 2007 but
% it would be fairly simple to extend to other datasets.
%
% As a first step, it is necessary to download the PASCAL VOC 2007 dataset.
% Go to http://pascallin.ecs.soton.ac.uk/challenges/VOC/voc2007/ and get:
%  VOCdevkit_08-Jun-2007, VOCtrainval_06-Nov-2007, VOCtest_06-Nov-2007
% After downloading these files extract them to dataDir='boxes/VOCdevkit'.
% Once complete this function will process the data and return information
% in a convenient format for boxesEval.m.
%
% USAGE
%  data = boxesData( opts )
%
% INPUTS
%  opts       - parameters (struct or name/value pairs)
%   .resDir     - ['boxes/'] location for results and evaluation
%   .dataDir    - ['boxes/VOCdevkit/'] dir containing PASCAL VOC 2007
%   .split      - ['val'] data split for evaluation
%
% OUTPUTS
%  data       - dataset split information
%   .split      - data split for evaluation
%   .n          - number of images
%   .ids        - list of string ids
%   .imgs       - list of image filenames
%   .gt         - ground truth for each image
%
% EXAMPLE
%
% See also edgeBoxesDemo, edgeBoxes, boxesEval
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar and Larry Zitnick, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

% get default parameters (unimportant parameters are undocumented)
dfs={ 'resDir','boxes/', 'dataDir','boxes/VOCdevkit/', 'split','val' };
o=getPrmDflt(varargin,dfs,1);

% locations of PASCAL VOC dataset
if(~exist(o.dataDir,'dir')), error('dataset not found, see help'); end
imDir=[o.dataDir 'VOC2007/JPEGImages/'];
gtDir=[o.dataDir 'VOC2007/Annotations/'];
imList=[o.dataDir 'VOC2007/ImageSets/Main/' o.split '.txt'];
addpath(genpath([o.dataDir 'VOCcode']));

% check if data already exists, if yes load and return
dataNm=[o.resDir '/GroundTruth' '-' o.split '.mat'];
if(exist(dataNm,'file')), data=load(dataNm); data=data.data; return; end

% get list of image ids
if(~exist(imList,'file')), error('ids file not found'); end
f=fopen(imList); ids=textscan(f,'%s %*s'); ids=ids{1}; fclose(f);

% generate list of image and gt filenames then load gt
n=length(ids); imgs=cell(n,1); gt=cell(n,1);
for i=1:n, imgs{i}=[imDir ids{i} '.jpg']; end
for i=1:n, gt{i}=[gtDir ids{i} '.xml']; end
if(~exist(imgs{1},'file')), error('images not found'); end
if(~exist(gt{1},'file')), error('annotations not found'); end
parfor i=1:n, [~,gt{i}]=bbGt('bbLoad',gt{i},'format',1); end

% create output structure and cache to disk
data=struct('split',o.split,'n',n,'ids',{ids},'imgs',{imgs},'gt',{gt});
if(~exist(o.resDir,'dir')), mkdir(resDir); end; save(dataNm,'data');

end
