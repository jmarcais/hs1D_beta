%----------------------------------------------------
% TEST: BILININTERP vs. CUDAINTERP2
%----------------------------------------------------
% Gerald Buchgraber, 2009
% gerald.buchgraber@student.tugraz.at
% 
% Build bilininterp and cudainterp2 mex files first
% (using build_mex.m)
% 
%----------------------------------------------------
% 
% Original test script from Alexander Huth
% provided with cudainterp2 implementation:
% http://www.mathworks.com/matlabcentral/fileexchange/20248
% 
%----------------------------------------------------
% 
% Benchmarks GPU vs CPU bilinear interpolation.
% runs 2-D bilinear interpolations on random data using
% CUDAINTERP2, BILININTERP and INTERP2.
%
% SPEEDUP contains the relative speed of INTERP2 and CUDAINTERP2 or
% BILININTERP.
% Thus if SPEEDUP(i,j) = 5, CUDAINTERP2 or BILININTERP
% completed interpolation of a random 300x300 matrices to
% NSIZE(i) x NSIZE(i) in 1/10 the time of INTERP2.
% 
%----------------------------------------------------

close all; clear all; clc;

addpath('./cudainterp2');

orig_size		= 300;
new_sizes		= 600:50:1800;
nrepeats		= 3;

gpu_results         = zeros(1,numel(new_sizes));
gpu_results_bilin   = zeros(1,numel(new_sizes));
cpu_results         = zeros(1,numel(new_sizes));

data = rand(orig_size);
for ns = 1:numel(new_sizes)
	disp(sprintf('New data size: %d elements\n',new_sizes(ns)^2));
	
	gpu_temp = zeros(1,nrepeats);
    gpu_temp_bilin = zeros(1,nrepeats);
	cpu_temp = zeros(1,nrepeats);
	
	for t = [1 1:nrepeats]
		tic;
    	temp = cudainterp2(data,new_sizes(ns),new_sizes(ns));
		gpu_temp(t) = toc;
        
        tic;
        temp = bilininterp(data,new_sizes(ns),new_sizes(ns));
		gpu_temp_bilin(t) = toc;
		
		tic;
		[x y] = meshgrid(linspace(1,orig_size,new_sizes(ns)),linspace(1,orig_size,new_sizes(ns)));
		temp = interp2(data,x,y,'*linear');
		cpu_temp(t) = toc;
	end
	
	gpu_results(ns) = mean(gpu_temp);
    gpu_results_bilin(ns) = mean(gpu_temp_bilin);
	cpu_results(ns) = mean(cpu_temp);
end

speedup_bilininterp = cpu_results./gpu_results_bilin;
speedup_cudainterp2 = cpu_results./gpu_results;
nsize = new_sizes.^2;

figure;
plot(nsize, speedup_bilininterp, 'r-*', nsize, speedup_cudainterp2, 'b-*');
title('Bilinear Interpolation (GPU vs CPU)');
xlabel('New data size (# of elements)');
ylabel('Relative speed (GPU/CPU)');
legend('Speedup: BILININTERP', 'Speedup: CUDAINTERP2', 'Location', 'NorthWest');
