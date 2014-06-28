% STRUCTUREDFOREST
% See also
%
% Fast edged detector code is based on the paper:
%  P. Dollár and C. Zitnick
%  "Structured Forests for Fast Edge Detection", ICCV 2013.
% Please cite the above paper if you end up using the edge detector.
% Code written and maintained by Piotr Dollar.
%
% Structured Edge detector code:
%   edgesChns     - Compute features for structured edge detection.
%   edgesDemo     - Demo for Structured Edge Detector (please see readme.txt first).
%   edgesDemoRgbd - Demo for RGBD Structured Edge Detector (please see readme.txt first).
%   edgesDetect   - Detect edges in image.
%   edgesSweeps   - Parameter sweeps for structured edge detector.
%   edgesTrain    - Train structured edge detector.
%
% Edge detection evaluation code:
%   edgesEval     - Run and evaluate structured edge detector on BSDS500.
%   edgesEvalDir  - Calculate edge precision/recall results for directory of edge images.
%   edgesEvalImg  - Calculate edge precision/recall results for single edge image.
%   edgesEvalPlot - Plot edge precision/recall results for directory of edge images.
