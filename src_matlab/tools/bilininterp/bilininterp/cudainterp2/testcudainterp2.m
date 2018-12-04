function [speedup nsize] = testcudainterp2()
%TESTCUDAINTERP2 Benchmarks GPU vs CPU bilinear interpolation.
%   [SPEEDUP DSIZE KSIZE] = TESTCUDACONV() runs 2-D bilinear
%   interpolations on random data using both CUDAINTERP2 and INTERP2.
%
%   SPEEDUP contains the relative speed of CONV2 and CUDACONV.
%   Thus if SPEEDUP(i,j) = 5, CUDAINTERP2 completed interpolation
%   of a random 300x300 matrices to NSIZE(i) x NSIZE(i) in 1/10 
%   the time of INTERP2.
%
%   NSIZE is the sqrt of the interpolated data size.

orig_size		= 300;
new_sizes		= 600:50:1800;
nrepeats		= 3;

gpu_results		= zeros(1,numel(new_sizes));
cpu_results		= zeros(1,numel(new_sizes));

data = rand(orig_size);
for ns = 1:numel(new_sizes)
	disp(sprintf('New data size: %d elements\n',new_sizes(ns)^2));
	
	gpu_temp = zeros(1,nrepeats);
	cpu_temp = zeros(1,nrepeats);
	
	for t = [1 1:nrepeats]
		tic;
		temp = cudainterp2(data,new_sizes(ns),new_sizes(ns));
		gpu_temp(t) = toc;
		
		tic;
		[x y] = meshgrid(linspace(1,orig_size,new_sizes(ns)),linspace(1,orig_size,new_sizes(ns)));
		temp = interp2(data,x,y,'*linear');
		cpu_temp(t) = toc;
	end
	
	gpu_results(ns) = mean(gpu_temp);
	cpu_results(ns) = mean(cpu_temp);
end

speedup = cpu_results./gpu_results;
nsize = new_sizes.^2;

figure;
plot(nsize,speedup,'-*');
xlabel('New data size (# of elements)')
ylabel('Relative speed (GPU/CPU)')

title('Bilinear Interpolation (GPU vs CPU)')

%dlmwrite('performance2.dat',speedup,',');