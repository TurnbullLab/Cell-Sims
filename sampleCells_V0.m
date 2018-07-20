% Set path to include binornd

% Seed random number generator using the clock
[seed0, seed0n] = resetRand;
% Can restart with same random seed:
% resetRand(seed0, seed0n);

debugPlot = true;  % Debug output (E.g. figure 1)

nCells = 10000;  % number of original ancestor cells to simulate
% sigma_t = 0.05;   % replace by sigma(t), prob of cell division in one time step (e.g. hour)
% alphaOverSig_slope = 0.01;  % replace by alpha(t)/sigma(t)   
HMax = 24 * 10;   % number of hours.

if debugPlot
  figure(1); clf; hold on;
  nPlot = 100;  % Only plot first nPlot samples 
  nPlot = min(nPlot, nCells);
end

mnN = zeros(1, HMax);
mnMhat = zeros(1, HMax);
stdN = zeros(1, HMax);
stdMhat = zeros(1, HMax);
N = ones(1, nCells);
Mhat = zeros(1, nCells);
for t = 1:HMax  
  alphaOverSig_t = 1-min(( * t, 1);
  idCell = N > 0;
  S = zeros(1, nCells);
  if any(idCell)
    S(idCell) = binornd(N(idCell), sigma_t);
  end
  A = zeros(1, nCells);
  idCell = S > 0;
  if any(idCell)
    A(idCell) = binornd(S(idCell), alphaOverSig_t);
  end
  N = N - S + 2*A;
  Mhat = Mhat + (S - A);
  
  figure(1); hold on;
  plot(repmat(t, 1, nPlot), Mhat(1:nPlot), '*b');
  plot(repmat(t, 1, nPlot), N(1:nPlot), 'og');

  mnN(t) = mean(N);
  mnMhat(t) = mean(Mhat);
  stdN(t) = std(N);
  stdMhat(t) = std(Mhat);

end

figure(2); clf;
plot(1:HMax, mnN, '-g');
hold on;
plot(1:HMax, mnMhat, '-b');
plot(1:HMax, mnN+stdN, '-c');
plot(1:HMax, mnN-stdN, '-c');

plot(1:HMax, mnMhat + stdMhat, '-m');
plot(1:HMax, mnMhat - stdMhat, '-m');
