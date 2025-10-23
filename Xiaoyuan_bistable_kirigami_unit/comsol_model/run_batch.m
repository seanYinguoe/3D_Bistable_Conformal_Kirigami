deltas = 0:0.1:8.5;

Ecurve = zeros(size(deltas));
theta_guess = [];  % warm start seed

for i = 1:numel(deltas)
    delta = deltas(i);

    % fast call: no plotting + provide warm start guess
    [~, Ei, theta_opt] = deform_triangle_batch(delta, edgeLen, l1, l4, beta, t, ...
        'theta_guess', theta_guess, 'do_plot', false);

    Ecurve(i)    = Ei;
    theta_guess  = theta_opt;  % chain the warm start
end

% quick plot
figure; plot(deltas, Ecurve, 'LineWidth',1.6);
xlabel('\delta'); ylabel('Energy');
grid on; title('Delta–Energy curve with warm start');
