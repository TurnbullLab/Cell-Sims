clear;              % clear parameters from other programs

debugPlot = true;   % Debug output (E.g. figure 1)

nCells = 1000;     % number of starting gcp cells to simulate

% Set up "par" structure for parameters needed in functions
par.Td      = 20;   % gcp doubling time (in h)
par.tau     = 8*24; % clone differentiation time (in h)
par.lambda  = 32;    % clone "tuning" constant (dimensionless)

% Time is specified in hours, from t0 to t1 with time step, dt
t0 = 2*24;          % start time(0=E15=P-4; 2 days = E17) 
t1 = 20*24;         % stop time (E15 + 20 days = P16)
dt = 1;             % time step (= 1h)

ts=t0:dt:t1;        % vector to hold times
kMax=length(ts);    % number of time steps

if debugPlot
  figure(1); clf; hold on;
  nPlot = 10;  % Only plot first nPlot samples 
  nPlot = min(nPlot, nCells);
end

mnN = zeros(1,kMax);        % mean of N 
mnMhat = zeros(1, kMax);    % mean of Mhat
stdN = zeros(1, kMax);      % stddev of N
stdMhat = zeros(1, kMax);   % stddev of Mhat

N = ones(1, nCells);        % N
Mhat = zeros(1, nCells);    % Mhat

for k=1:kMax
    t=ts(k);                % time is now indexed by k
    idCell = N > 0;
    S = zeros(1, nCells);
    if any(idCell)
        S(idCell) = binornd(N(idCell), sigma(t,par)*dt);
    end
    A = zeros(1, nCells);
    idCell = S > 0;
    if any(idCell)
        A(idCell) = binornd(S(idCell), alpha(t,par)/sigma(t,par));
    end
    
  N = N - S + 2*A;
  Mhat = Mhat + (S - A);
  
  figure(1); hold on;       % Plots sample simulated gcp clones
  plot(repmat(t, 1, nPlot), 2*Mhat(1:nPlot), '.b'); % #gcs = 2*Mhat
  plot(repmat(t, 1, nPlot), N(1:nPlot), '.g');      % #gcps= N

  mnN(k) = mean(N);         % Mean #gcps
  mnMhat(k) = 2*mean(Mhat); % Mean #gcs = 2*mean(Mhat)
  stdN(k) = std(N);         % Stdev #gcps
  stdMhat(k) = 2*std(Mhat); % Stdev #gcs = 2*std(Mhat)
 
 % Use this to examine histograms of N and MHat 
 % if mod(k, 200)==0
 %    figure(3); clf;
 %    histogram(N, 100, 'FaceColor', 'g');
 %    hold on;
 %    h = histogram(2*Mhat, 100, 'FaceColor', 'b');
 %    pause;
 % end
end

hCurve = [];
fontSize = 16;
lineWidth = 3;
figure(2); clf;             % Plots mean+/-stdev #'s of gcps and gcs 
set(gca, 'FontSize', fontSize);
hCurve(1) = plot(ts(1:kMax), mnN, '-g', 'LineWidth', lineWidth);
hold on;
hCurve(2) = plot(ts(1:kMax), mnMhat, '-b', 'LineWidth', lineWidth);
hCurve(3) = plot(ts(1:kMax), mnN+stdN, '-c', 'LineWidth', lineWidth-1);
plot(ts(1:kMax), mnN-stdN, '-c', 'LineWidth', lineWidth-1);

hCurve(4) = plot(ts(1:kMax), mnMhat + stdMhat, '-m', 'LineWidth', lineWidth-1);
plot(ts(1:kMax), mnMhat - stdMhat, '-m' , 'LineWidth', lineWidth-1);

hl = legend(hCurve, 'N', 'M', 'N \pm std', 'M \pm std', 'Location', ...
       'NorthWest');
set(hl, 'FontSize', fontSize);
title('Cell Counts');
xlabel('Time');
ylabel('Number of Cells');

% You can then resize the figure by hand and export it (using
% File->Export Setup -> Export and choosing a postscript file
% output.  (It doesn't give identical results).

% Or you can try to print this as a postscript file 'foo.eps' directly:
%print(figure(2), 'depsc', 'foo');  
% But using print like this is far from  WYSIWYG, which is
% annoying.

% Or resize the figure by hand and take a screenshot....
%plot(ts(1:kMax), mnN, '-g');
%hold on;
%plot(ts(1:kMax), mnMhat, '-b');

%plot(ts(1:kMax), mnN+stdN, '-c');
%plot(ts(1:kMax), mnN-stdN, '-c');

%plot(ts(1:kMax), mnMhat+stdMhat, '-m');
%plot(ts(1:kMax), mnMhat-stdMhat, '-m');
