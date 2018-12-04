%----------------------------------------------------
% Gerald Buchgraber, 2009
% gerald.buchgraber@student.tugraz.at
%----------------------------------------------------
% NVCC MEX BUILD SCRIPT
% with updated options file for Matlab r2008a from:
% http://forums.nvidia.com/index.php?showtopic=66961
%
% PLEASE ENSURE VALID INCLUDE AND LIB PATHS BELOW!
%----------------------------------------------------

clear all; clc;

% compiling BILININTERP
nvmex -f nvmexopts_r2008a.bat bilininterp.cu -IC:\cuda\include -LC:\cuda\lib -lcudart

% compiling CUDAINTERP2
nvmex -f nvmexopts_r2008a.bat ./cudainterp2/cudainterp2.cu -IC:\cuda\include -I"C:\Program Files\NVIDIA Corporation\NVIDIA CUDA SDK\common\inc" -LC:\cuda\lib -lcudart -outdir ./cudainterp2