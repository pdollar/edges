###################################################################
#                                                                 #
#    Structured Edge Detection Toolbox V3.0                       #
#    Piotr Dollar (pdollar-at-gmail.com)                          #
#                                                                 #
###################################################################

1. Introduction.

Very fast edge detector (up to 60 fps depending on parameter settings) that achieves excellent accuracy. Can serve as input to any vision algorithm requiring high quality edge maps. Toolbox also includes the Edge Boxes object proposal generation method and fast superpixel code.

If you use the Structured Edge Detection Toolbox, we appreciate it if you cite an appropriate subset of the following papers:

@inproceedings{DollarICCV13edges,
  author    = {Piotr Doll\'ar and C. Lawrence Zitnick},
  title     = {Structured Forests for Fast Edge Detection},
  booktitle = {ICCV},
  year      = {2013},
}

@article{DollarARXIV14edges,
  author    = {Piotr Doll\'ar and C. Lawrence Zitnick},
  title     = {Fast Edge Detection Using Structured Forests},
  journal   = {ArXiv},
  year      = {2014},
}

@inproceedings{ZitnickECCV14edgeBoxes,
  author    = {C. Lawrence Zitnick and Piotr Doll\'ar},
  title     = {Edge Boxes: Locating Object Proposals from Edges},
  booktitle = {ECCV},
  year      = {2014},
}

###################################################################

2. License.

This code is published under the MSR-LA Full Rights License.
Please read license.txt for more info.

###################################################################

3. Installation.

a) This code is written for the Matlab interpreter (tested with versions R2013a-2013b) and requires the Matlab Image Processing Toolbox. 

b) Additionally, Piotr's Matlab Toolbox (version 3.26 or later) is also required. It can be downloaded at:
 https://pdollar.github.io/toolbox/.

c) Next, please compile mex code from within Matlab (note: win64/linux64 binaries included):
  mex private/edgesDetectMex.cpp -outdir private [OMPPARAMS]
  mex private/edgesNmsMex.cpp    -outdir private [OMPPARAMS]
  mex private/spDetectMex.cpp    -outdir private [OMPPARAMS]
  mex private/edgeBoxesMex.cpp   -outdir private
Here [OMPPARAMS] are parameters for OpenMP and are OS and compiler dependent.
  Windows:  [OMPPARAMS] = '-DUSEOMP' 'OPTIMFLAGS="$OPTIMFLAGS' '/openmp"'
  Linux V1: [OMPPARAMS] = '-DUSEOMP' CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
  Linux V2: [OMPPARAMS] = '-DUSEOMP' CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
To compile without OpenMP simply omit [OMPPARAMS]; note that code will be single threaded in this case.

d) Add edge detection code to Matlab path (change to current directory first): 
 >> addpath(pwd); savepath;

e) Finally, optionally download the BSDS500 dataset (necessary for training/evaluation):
 http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/
 After downloading BSR/ should contain BSDS500, bench, and documentation.

f) A fully trained edge model for RGB images is available as part of this release. Additional models are available online, including RGBD/D/RGB models trained on the NYU depth dataset and a larger more accurate BSDS model.

###################################################################

4. Getting Started.

 - Make sure to carefully follow the installation instructions above.
 - Please see "edgesDemo.m", "edgeBoxesDemo" and "spDemo.m" to run demos and get basic usage information.
 - For a detailed list of functionality see "Contents.m".

###################################################################

5. History.

Version NEW
 - now hosting on github (https://github.com/pdollar/edges)
 - suppress Mac warnings, added Mac binaries
 - edgeBoxes: added adaptive nms variant described in arXiv15 paper

Version 3.01 (09/08/2014)
 - spAffinities: minor fix (memory initialization)
 - edgesDetect: minor fix (multiscale / multiple output case)

Version 3.0 (07/23/2014)
 - added Edge Boxes code corresponding to ECCV paper
 - added Sticky Superpixels code
 - edge detection code unchanged

Version 2.0 (06/20/2014)
 - second version corresponding to arXiv paper
 - added sharpening option
 - added evaluation and visualization code
 - added NYUD demo and sweep support
 - various tweaks/improvements/optimizations

Version 1.0 (11/12/2013)
 - initial version corresponding to ICCV paper

###################################################################
