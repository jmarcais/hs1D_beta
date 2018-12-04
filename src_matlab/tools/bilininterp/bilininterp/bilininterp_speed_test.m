%----------------------------------------------------
% SPEED TEST
%----------------------------------------------------
% Gerald Buchgraber, 2009
% gerald.buchgraber@student.tugraz.at
% 
% Build bilininterp mex file first
% (using build_mex.m)
%----------------------------------------------------

close all; clear all; clc;

%----------------------------------------------------
% C O N F I G
%----------------------------------------------------
basic_grid_elements_x = 256;
basic_grid_elements_y = 256;

test_trajectories = round(linspace(1, 10000, 10));
samples_on_each_trajectory = 300;
num_avg_repeats = 5;
%----------------------------------------------------

bilininterp_results = zeros(1, numel(test_trajectories));
interp2_results = zeros(1, numel(test_trajectories));

for i = 1:numel(test_trajectories)
    numSpokes = test_trajectories(i);
    numSamplesOnSpoke = samples_on_each_trajectory;

    theta = 0:pi/numSpokes:pi-pi/numSpokes;
    rho = linspace(-0.5, 0.5, numSamplesOnSpoke)';

    b = zeros(numSamplesOnSpoke, numSpokes);
    for ii = 1:length(theta)
        b(:,ii) = rho*exp(-j*theta(ii));
    end
    clear theta;

    XI = real(b);
    YI = imag(b);

    fx = round(basic_grid_elements_x / 3);
    fy = round(basic_grid_elements_y / 3);
    [X, Y] = meshgrid(1:1:3*fx, 1:1:3*fy);
    Z = peaks(2*X/fx-3, 2*Y/fy-3);

    XI = (XI + .5) * (max(max(X))-1) + 1;
    YI = (YI + .5) * (max(max(Y))-1) + 1;
    
    XI = sort(XI(:));
    YI = sort(YI(:));
    
    bilininterp_temp = zeros(1, num_avg_repeats);
    interp2_temp = zeros(1, num_avg_repeats);
    for t = 1:num_avg_repeats
        %----------------------------------------------------
        % BILININTERP
        %----------------------------------------------------
        tic;
        ZBI = bilininterp(Z, XI, YI);
        bilininterp_temp(t) = toc;
        clear ZBI;
    end
    for t = 1:num_avg_repeats
        %----------------------------------------------------
        % INTERP2
        %----------------------------------------------------
        tic;
        ZII = interp2(Z, XI, YI, '*linear');
		interp2_temp(t) = toc;
        clear ZII;
    end
    bilininterp_results(i) = mean(bilininterp_temp);
    interp2_results(i) = mean(interp2_temp);
    clc; disp(['Progress: ', num2str(i / numel(test_trajectories) * 100), '%']);
end

% Visualize results
%----------------------------------------------------
figure;
samples = test_trajectories * samples_on_each_trajectory;
plot(samples, bilininterp_results, 'r-*', samples, interp2_results, 'g-*');
title('Speed Test Results');
xlabel('# of interpolated samples');
ylabel('Time [s]');
legend('Bilinear Interpolation based on CUDA', 'Matlab function: interp2', 'Location', 'NorthWest');

% Speedup
%----------------------------------------------------
speedup = interp2_results ./ bilininterp_results;
nsize = samples;
figure;
plot(nsize, speedup, '-*');
title('SPEEDUP - bilininterp vs. interp2');
xlabel('# of interpolated samples');
ylabel('Speedup: t_{CPU} / t_{GPU}');

