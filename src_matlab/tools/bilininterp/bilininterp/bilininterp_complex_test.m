%----------------------------------------------------
% TEST: BILINEAR INTERPOLATION WITH COMPLEX NUMBERS
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
orgGridNumRows = 10;
orgGridNumCols = 10;

[X,Y] = meshgrid(1:orgGridNumCols, 1:orgGridNumRows);
Z = complex(rand(orgGridNumRows, orgGridNumCols), rand(orgGridNumRows, orgGridNumCols)) * 20;

interpGridNumRows = 30;
interpGridNumCols = 30;
XYI = ones(interpGridNumRows, interpGridNumCols);
stepsX = linspace(1, orgGridNumCols, interpGridNumCols);
stepsY = linspace(1, orgGridNumRows, interpGridNumRows);
for i = 1:numel(stepsX)
    for j = 1:numel(stepsY)
        XYI(j,i) = complex(stepsX(i), stepsY(j));
    end
end

%----------------------------------------------------
% BILININTERP
%----------------------------------------------------
ZBI = bilininterp(Z, real(XYI), imag(XYI));

%----------------------------------------------------
% INTERP2
%----------------------------------------------------
ZII = interp2(Z, real(XYI), imag(XYI), '*linear');


% Visualization
%----------------------------------------------------
figure;
mesh(real(XYI), imag(XYI), real(ZBI) + 40), hold;
mesh(real(XYI), imag(XYI), real(ZII) + 20);
mesh(X, Y, real(Z));
hold off;
title('Real part');
legend('CUDA: bilininterp', 'Matlab: interp2', 'Source', 'Location', 'NorthWest');


figure;
mesh(real(XYI), imag(XYI), imag(ZBI) + 40), hold;
mesh(real(XYI), imag(XYI), imag(ZII) + 20);
mesh(X, Y, imag(Z));
hold off;
title('Imaginary part');
legend('CUDA: bilininterp', 'Matlab: interp2', 'Source', 'Location', 'NorthWest');


figure;
mesh(real(XYI), imag(XYI), abs(ZBI - ZII)./abs(ZII));
title('Relative error');
