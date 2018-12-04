%BILININTERP 2-D bilinear interpolation using graphic hardware
%   B = BILININTERP(A, NR, NC) performs a bilinear interpolation
%   on the matrix A so that the output B has NR rows and NC cols.
%
%   The interpolation is based on Nvidia CUDA.
%   
%   Class support for input A: 
%      float (real): double
%
%   Class support for inputs NR, NC: 
%      real
    
%   By Gerald Buchgraber, Technical University of Graz, Austria, 2009
%   inspired by CUDAINTERP2 (Alexander Huth)
%   http://www.mathworks.com/matlabcentral/fileexchange/20248
