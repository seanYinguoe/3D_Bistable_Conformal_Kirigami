% parfor_vs_for_demo.m
% Run this script from MATLAB. Requires Parallel Computing Toolbox.

clear; clc;

N = 2e5;                         % increase for heavier test
x = rand(N,1);

% Test function: CPU-heavy independent work per iteration
work = @(v) sum(sin(v*(1:200)).^2 + cos(v*(1:200)).^2);

% --- Serial for ---
y_for = zeros(N,1);
tic;
for i = 1:N
    y_for(i) = work(x(i));
end
t_for = toc;

% --- Parallel parfor ---
y_par = zeros(N,1);
tic;
parfor i = 1:N
    y_par(i) = work(x(i));
end
t_par = toc;

% Check correctness and speedup
max_err = max(abs(y_for - y_par));
speedup = t_for / t_par;

fprintf('for time:    %.3f s\n', t_for);
fprintf('parfor time: %.3f s\n', t_par);
fprintf('speedup:     %.2fx\n', speedup);
fprintf('max error:   %.3e\n', max_err);
