function TestScaleTime(doSpeed)
% Test ScaleTime
%
% Tested: Matlab 6.5, 7.7
% Author: Jan Simon, Heidelberg, (C) 2009 J@n-Simon.De

% $JRev: R0f V:005 Sum:265F0599 Date:30-Sep-2009 01:00:47 $
% $File: User\JSim\Published\ScaleTime\TestScaleTime.m $

% Initialize: ==================================================================
ErrID = ['JSim:', mfilename];
whichScaleTime = which('ScaleTime');

disp(['==== Test ScaleTime:  ', datestr(now, 0), char(10), ...
      'Version: ', whichScaleTime, char(10)]);

% Key for linear interpolation - the source of INTERP1 in Matlab6.5 states, that
% "*linear" is the old Matlab 5.3 method. But for large arrays the test of
% equally spaced steps is time-consuming!
% Set it to 'linear', if '*linear' does not work anymore.
linearMethod = '*linear';

if nargin == 0
   doSpeed = true;
end

% Look for lininterp1f, Umberto Picchini, Nov 2005
% http://www.mathworks.com/matlabcentral/fileexchange/8627
hasInterp1f = ~isempty(which('lininterp1f'));

% Look for qinterp1, Nathaniel Brahms,
% http://www.mathworks.com/matlabcentral/fileexchange/10286
hasQinterp1 = ~isempty(which('qinterp1'));

% Do we use the MEX or M version:
[dummy1, dummy2, fileExt] = fileparts(whichScaleTime);  %#ok<ASGLU>
useMex = strcmpi(strrep(fileExt, '.', ''), mexext);

% Start tests: -----------------------------------------------------------------
try
   desc = 'ScaleTime([], [])';
   y = ScaleTime([], []);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isempty(y)
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime([], 1, 2, 0)';
   y = ScaleTime([], 1, 2, 0);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isempty(y)
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime([], 1, 2, -1)';
   y = ScaleTime([], 1, 2, -1);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isempty(y)
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime([1; 2], 1.5, 1, 1)';
   y = ScaleTime([1; 2], 1.5, 1, 1);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isempty(y)
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime((1:10)'', [])';
   y = ScaleTime((1:10)', []);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isempty(y)
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime(ones(10, 10), [])';
   y = ScaleTime(ones(10, 10), []);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isempty(y)
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime((1:10)'', 1)';
   y = ScaleTime((1:10)', 1);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isequal(y, 1)
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime(ones(10, 10), 1)';
   y = ScaleTime(ones(10, 10), 1);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isequal(y, ones(1, 10))
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime(ones(10, 10), 10)';
   y = ScaleTime(ones(10, 10), 10);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isequal(y, ones(1, 10))
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime([1; 2], 1, 1, 1)';
   y = ScaleTime([1; 2], 1, 1, 1);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isequal(y, 1)
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime([1; 2], 1.5, 2, 1)';
   y = ScaleTime([1; 2], 1.5, 2, 1);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, 1.5)
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime([1; 3], 1.6)';
   y = ScaleTime([1; 3], 1.6);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, 2.2)
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime((1:10)'', 1, 10, 9)';
   y = ScaleTime((1:10)', 1, 10, 9);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, (1:1.125:10)')
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime((1:10)'', 1, 10, 10)';
   y = ScaleTime((1:10)', 1, 10, 10);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, (1:10)')
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime((1:10)'', 1, 10, 11)';
   y = ScaleTime((1:10)', 1, 10, 11);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, (1:0.9:10)')
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime((1:10)'', 1, 9, 8)';
   y = ScaleTime((1:10)', 1, 9, 8);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, (1:(8 / 7):9)')
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime((1:10)'', 1, 9, 9)';
   y = ScaleTime((1:10)', 1, 9, 9);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, (1:9)')
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime((1:10)'', 1, 9, 10)';
   y = ScaleTime((1:10)', 1, 9, 10);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, (1:8/9:9)')
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime(1:10, 1, 9, 10) row vector';
   data = rand(1, 10);
   yRow = ScaleTime(data,    1, 9, 10);
   yCol = ScaleTime(data(:), 1, 9, 10);
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isequal(yRow(:), yCol(:)) && isequal(size(yRow), [1, 10])
   disp(['  ok: ', desc]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

% Compare with INTERP1: --------------------------------------------------------
disp([char(10), '== Compare with INTERP1:']);
try
   desc = 'ScaleTime(sin, 1, 2*pi, 100)';
   data = sin(0:0.1:2*pi)';
   y    = ScaleTime(data, 1, 2*pi, 100);
   yM   = interp1(1:size(data, 1), data, linspace(1, 2*pi, 100)');
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, yM)
   disp(['  ok: ', desc, ' == INTERP1']);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime([sin, sin(down:up)], 1, 2*pi, 100)';
   data = cat(2, data, data(end:-1:1));
   y    = ScaleTime(data, 1, 2*pi, 100);
   yM   = interp1(1:size(data, 1), data, linspace(1, 2*pi, 100)');
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, yM)
   disp(['  ok: ', desc, ' == INTERP1']);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

try
   desc = 'ScaleTime(rand(100, 100), 1, 100, 57)';
   data = rand(100, 100);
   y    = ScaleTime(data, 1, 100, 57);
   yV   = ScaleTime(data, 1:99/56:100);
   yM   = interp1(1:size(data, 1), data, linspace(1, 100, 57)');
catch
   error(ErrID, ['Crash: ', desc, char(10), lasterr]);
end
if isEqualTol(y, yM, 100 * eps)
   disp(['  ok: ', desc, ' == INTERP1']);
   disp(['      max deviation: ', ...
         sprintf('%.1f * eps', max(abs(y(:) - yM(:))) / eps)]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end
if isEqualTol(y, yV, 100 * eps)
   disp(['  ok: ', desc, ' == ScaleTime(rand(100, 100), 1:99/56:100)']);
   disp(['      max deviation: ', ...
         sprintf('%.1f * eps', max(abs(yM(:) - yV(:))) / eps)]);
else
   error(ErrID, ['Bad reply from: ', desc]);
end

% Speed: -----------------------------------------------------------------------
if doSpeed
   TestTime = 0.5;  % sec
   fprintf('\n== Speed test (test time: %g sec):\n', TestTime);
else
   TestTime = 0.1;
   fprintf('\n== Speed test (test time: %g sec - may be inaccurate):\n', ...
      TestTime);
end

DataLen   = [10, 100, 1000, 10000, 100000];
DataWidth = [1, 10, 100];
InterpLen = [10, 100, 1000, 10000];

for aDataLen = DataLen
   for aDataWidth = DataWidth
      fprintf('Data size: [%d x %d]\n', aDataLen, aDataWidth);
      fprintf('  Interpolation steps:  ');
      fprintf('  %7d', InterpLen);
      fprintf('\n    INTERP1:            ');
      for aInterpLen = InterpLen
         Data = rand(aDataLen, aDataWidth);
         NewT = 1:(aDataLen - 1) / (aInterpLen - 1):aDataLen;
         
         % At least 2 loops!
         iTime     = cputime;
         y         = interp1(1:aDataLen, Data, NewT, linearMethod); %#ok<*NASGU>
         y         = interp1(1:aDataLen, Data, NewT, linearMethod);
         interp1_N = 2;
         while cputime - iTime < TestTime
            y         = interp1(1:aDataLen, Data, NewT, linearMethod);
            interp1_N = interp1_N + 1;
         end
         aTime        = cputime - iTime;
         interp1_Loop = interp1_N / aTime;  % Loops per second
         PrintLoop(interp1_Loop);
      end
      fprintf('  loops/sec\n');
      drawnow;
      
      fprintf('    ScaleTime (vector): ');
      for aInterpLen = InterpLen
         Data = rand(aDataLen, aDataWidth);
         NewT = 1:(aDataLen - 1) / (aInterpLen - 1):aDataLen;
         
         iTime    = cputime;
         y        = ScaleTime(Data, NewT);
         y        = ScaleTime(Data, NewT);
         ScaleV_N = 2;
         while cputime - iTime < TestTime
            y        = ScaleTime(Data, NewT);
            ScaleV_N = ScaleV_N + 1;
         end
         aTime       = cputime - iTime;
         ScaleV_Loop = ScaleV_N / aTime;  % Loops per second
         PrintLoop(ScaleV_Loop);
      end
      fprintf('  loops/sec\n');
      drawnow;
      
      fprintf('    ScaleTime (index):  ');
      for aInterpLen = InterpLen
         iTime    = cputime;
         y        = ScaleTime(Data, 1, aDataLen, aInterpLen);
         y        = ScaleTime(Data, 1, aDataLen, aInterpLen);
         ScaleT_N = 2;
         while cputime - iTime < TestTime
            y        = ScaleTime(Data, 1, aDataLen, aInterpLen);
            ScaleT_N = ScaleT_N + 1;
         end
         aTime       = cputime - iTime;
         ScaleT_Loop = ScaleT_N / aTime;  % Loops per second
         PrintLoop(ScaleT_Loop);
      end
      fprintf('  loops/sec\n');
      drawnow;
      
      % Compare with local Matlab version:
      if useMex
         fprintf('    ScaleTime (Matlab): ');
         for aInterpLen = InterpLen
            iTime    = cputime;
            NewT     = 1:(aDataLen - 1) / (aInterpLen - 1):aDataLen;
            y        = ScaleTime_local(Data, NewT);
            y        = ScaleTime_local(Data, NewT);
            ScaleM_N = 2;
            while cputime - iTime < TestTime
               y        = ScaleTime_local(Data, NewT);
               ScaleM_N = ScaleM_N + 1;
            end
            aTime       = cputime - iTime;
            ScaleM_Loop = ScaleM_N / aTime;  % Loops per second
            PrintLoop(ScaleM_Loop);
         end
         fprintf('  loops/sec\n');
         drawnow;
      end
      
      % Compare with functions from Matlab's FEX:
      if hasInterp1f && aDataWidth == 1
         fprintf('    lininterp1f:        ');
         for aInterpLen = InterpLen
            iTime      = cputime;
            NewT       = 1:(aDataLen - 1) / (aInterpLen - 1):aDataLen;
            y          = lininterp1f(1:length(Data), Data, NewT, []);
            y          = lininterp1f(1:length(Data), Data, NewT, []);
            interp1f_N = 2;
            while cputime - iTime < TestTime
               y          = lininterp1f(1:length(Data), Data, NewT, []);
               interp1f_N = interp1f_N + 1;
            end
            aTime         = cputime - iTime;
            interp1f_Loop = interp1f_N / aTime;  % Loops per second
            PrintLoop(interp1f_Loop);
         end
         fprintf('  loops/sec\n');
         drawnow;
      end
      
      if hasQinterp1 && aDataWidth == 1
         fprintf('    qinterp1:           ');
         for aInterpLen = InterpLen
            iTime      = cputime;
            NewT       = 1:(aDataLen - 1) / (aInterpLen - 1):aDataLen;
            y          = qinterp1(1:length(Data), Data, NewT, 1);
            y          = qinterp1(1:length(Data), Data, NewT, 1);
            qinterp1_N = 2;
            while cputime - iTime < TestTime
               y          = qinterp1(1:length(Data), Data, NewT, 1);
               qinterp1_N = qinterp1_N + 1;
            end
            aTime         = cputime - iTime;
            qinterp1_Loop = qinterp1_N / aTime;  % Loops per second
            PrintLoop(qinterp1_Loop);
         end
         fprintf('  loops/sec\n');
         drawnow;
      end
      
      fprintf('\n');
   end
end

return;

% ******************************************************************************
function PrintLoop(N)
if N > 10
   fprintf('  %7.0f', N);
else
   fprintf('  %7.1f', N);
end

return;

% ******************************************************************************
% Local copy of ScaleTime to test speed if compiled MEX is in the path:
function Yi = ScaleTime_local(Y, Ti, Tf, Tn)
% Author: Jan Simon, Heidelberg, (C) 2009 J@n-Simon.De
% Was JRev: R0h V:033 Sum:98748BA6 Date:30-Sep-2009 00:50:13 $
% Was File: Tools\GLMath\ScaleTime.m $

[nRow, nCol] = size(Y);
doTranspose  = (nRow == 1 && nCol > 1);
if doTranspose
   nRow = nCol;
   nCol = 1;
   Y    = Y(:);
end

switch nargin
   case 4
      if and(Ti < 1, Tf > nRow)
         error(['*** ', mfilename, ': Interpolation frames out of range.']);
      end
      
      if Tn < 1.0
         Yi = [];
         return;
      elseif Tn > 1
         Ti = Ti:((Tf - Ti) / (Tn - 1)):Tf;
      elseif Tf < Ti
         Ti = [];
      end
      
   case 2
      if any(Ti < 1)
         error(['*** ', mfilename, ': All [Ti] must be >= 1.']);
      end
      if any(Ti > nRow)
         error(['*** ', mfilename, ...
               ': All [Ti] must be <= number of rows of Y.']);
      end
      
   otherwise
      error(['*** ', mfilename, ': 2 or 4 inputs required.']);
end

if isempty(Ti)
   Yi = [];
   return;
end

if nRow < 2
   error(['*** ', mfilename, ': Y must have at least 2 rows.']);
end

Ti = Ti(:);
Si = Ti - floor(Ti);
Ti = floor(Ti);

d     = (Ti == nRow);
Ti(d) = Ti(d) - 1;
Si(d) = 1;

if nCol > 1
   Si = Si(:, ones(1, nCol));
end
Yi = Y(Ti, :) .* (1 - Si) + Y(Ti + 1, :) .* Si;

if doTranspose
   Yi = reshape(Yi, 1, []);
end

return;
