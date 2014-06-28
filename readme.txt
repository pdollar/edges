###################################################################
#                                                                 #
#    Structured Edge Detection Toolbox V2.0                       #
#    Piotr Dollar (pdollar-at-microsoft.com)                      #
#                                                                 #
###################################################################

1. Introduction.

Very fast edge detector (1-60 fps depending on parameter settings) that achieves excellent accuracy (top accuracy on BSDS500 Segmentation dataset and NYU Depth dataset as of publication date). Can serve as input to any vision algorithm requiring high quality edge maps. 

If you use the Structured Edge Detection Toolbox, we appreciate it if you cite the following work in any resulting publication:

@inproceedings{DollarICCV13edges,
  author={Piotr Doll\'ar and C. Lawrence Zitnick},
  title={Structured Forests for Fast Edge Detection},
  booktitle={ICCV},
  year={2013},
}

@article{2014arXiv,
  author={Piotr Doll\'ar and C. Lawrence Zitnick},
  title={Fast Edge Detection Using Structured Forests},
  journal = {ArXiv},
  year = 2014,
}

###################################################################

2. License.

This code is published under the MSR-LA Full Rights License.
Please read license.txt for more info.

###################################################################

3. Installation.

a) This code is written for the Matlab interpreter (tested with versions R2013a-2013b) and requires the Matlab Image Processing Toolbox. 

b) Additionally, Piotr's Matlab Toolbox (version 3.24 or later) is also required. It can be downloaded at: 
 http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html.

c) Next, please compile mex code from within Matlab (note: win64/linux64 binaries included):
Windows:
  mex private/edgesDetectMex.cpp -outdir private '-DUSEOMP' 'OPTIMFLAGS="$OPTIMFLAGS' '/openmp"'
  mex private/edgesNmsMex.cpp    -outdir private '-DUSEOMP' 'OPTIMFLAGS="$OPTIMFLAGS' '/openmp"'
Linux version 1:
  mex private/edgesDetectMex.cpp -outdir private '-DUSEOMP' CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
  mex private/edgesNmsMex.cpp    -outdir private '-DUSEOMP' CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
Linux version 2:
  mex private/edgesDetectMex.cpp -outdir private '-DUSEOMP' CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
  mex private/edgesNmsMex.cpp    -outdir private '-DUSEOMP' CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

d) Add edge detection code to Matlab path (change to current directory first): 
 >> addpath(pwd); savepath;

e) Finally, optionally download the BSDS500 dataset (necessary for training/evaluation):
 http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/
 After downloading BSR/ should contain BSDS500, bench, and documentation.

f) A fully trained edge model for RGB images is available as part of this release. Additional models are available online, including RGBD/D/RGB models trained on the NYU depth dataset and a larger more accurate BSDS model.

###################################################################

4. Getting Started.

 - Make sure to carefully follow the installation instructions above.
 - Please see "edgesDemo.m" to run a demo and get basic usage information.
 - For a detailed list of functionality see "Contents.m".

###################################################################

5. History.

Version 2.0 (06/20/2014)
 - second version corresponding to arXiv paper
 - added sharpening option
 - added evaluation and visualization code
 - added NYUD demo and sweep support
 - various tweaks/improvements/optimizations

Version 1.0 (11/12/2013)
 - initial version corresponding to ICCV paper

###################################################################
