%----------------------------------------------------
% TEST: TRAJECTORY INTERPOLATION AND VISUALIZATION
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
numSpokes = 80;
numSamplesOnSpoke = 50;
% numSpokes = 500;
% numSamplesOnSpoke = 300;
%----------------------------------------------------

theta = 0:pi/numSpokes:pi-pi/numSpokes;
rho = linspace(-0.5, 0.5, numSamplesOnSpoke)';
b = zeros(numSamplesOnSpoke, numSpokes);
for ii = 1:length(theta)
    b(:,ii) = rho*exp(-j*theta(ii));
end
clear theta;

XI = real(b);
YI = imag(b);

[X, Y] = meshgrid(1:1:30);
Z = peaks(X/5-3, Y/5-3);

XI = (XI + .5) * (max(max(X))-1) + 1;
YI = (YI + .5) * (max(max(Y))-1) + 1;

% show boundary difference to interp2
% XI = XI - 10;
% YI = YI - 10;

%----------------------------------------------------
% BILININTERP
%----------------------------------------------------
tic;
ZBI = bilininterp(Z, XI, YI);
disp(['CUDA bilininterp: ', num2str(toc)]);

%----------------------------------------------------
% INTERP2
%----------------------------------------------------
tic;
ZII = interp2(Z, XI, YI, '*linear');
disp(['Matlab interp2: ', num2str(toc)]);

mesh(XI, YI, ZBI + 20), hold;
mesh(XI, YI, ZII + 10);
mesh(X, Y, Z);
hold off;
title('BILININTERP - Bilinear Interpolation Test');
legend('CUDA: bilininterp', 'Matlab: interp2', 'Source', 'Location', 'NorthWest');
