%% Parameters (Modify as needed)
k = 1;      % Spring stiffness
edgeLen = 15;
l1 = edgeLen * 0.8;
l2 = edgeLen * 0.1;
t  = edgeLen * 0.02;
l3 = l1 -2*l2 - (edgeLen - l1 -l2)/2;  % stretch length
delta_min = 0;  % Min displacement (>0)
delta_max = 10;  % Max displacement
num_points = 1000; % Resolution

%% Derived parameters
rest_alpha1 = pi/3;       % Rest angle α₁
rest_alpha2 = 2*pi/3;     % Rest angle α₂

%% Define functions
% Potential energy
U = @(a1, a2) 0.5*k*(a1 - rest_alpha1).^2 + 0.5*k*(a2 - rest_alpha2).^2;

% Constraints (f1=0, f2=-δ)
f1 = @(a1, a2) l2*(-sqrt(3) + sin(2*pi/3 - a1) + sin(a1)) + ...
                l3*(-sqrt(3)/2 + sin(a1 + a2 - 2*pi/3));
            
f2 = @(a1, a2) l2*(-1 - cos(2*pi/3 - a1) + cos(a1)) + ...
                l3*(-0.5 + cos(a1 + a2 - 2*pi/3)) + l2;

%% Part 1: Energy-Displacement Curve (U vs δ)
delta_vec = linspace(delta_min, delta_max, num_points);
alpha1_vec = zeros(size(delta_vec));
alpha2_vec = zeros(size(delta_vec));
U_vec = zeros(size(delta_vec));

% Initial guess (rest position)
a1_guess = rest_alpha1;
a2_guess = rest_alpha2;

options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-9);

for i = 1:length(delta_vec)
    delta_val = delta_vec(i);
    
    % Solve constraints: [f1=0, f2=-δ]
    constraint_eq = @(x) [f1(x(1), x(2)); 
                           f2(x(1), x(2)) + delta_val];
    
    [sol, ~, exitflag] = fsolve(constraint_eq, [a1_guess; a2_guess], options);
    
    if exitflag > 0
        alpha1_vec(i) = sol(1);
        alpha2_vec(i) = sol(2);
        U_vec(i) = U(sol(1), sol(2));
        a1_guess = sol(1); % Update guess for next iteration
        a2_guess = sol(2);
    else
        alpha1_vec(i) = NaN;
        alpha2_vec(i) = NaN;
        U_vec(i) = NaN;
    end
end

%% Part 2: Find Local Extrema (where λ₂=0)
% System: [f1=0, ∂U/∂α₁ - λ₁*∂f1/∂α₁=0, ∂U/∂α₂ - λ₁*∂f1/∂α₂=0]
df1_da1 = @(a1, a2) l2*(cos(a1) - cos(2*pi/3 - a1)) + ...
                     l3*cos(a1 + a2 - 2*pi/3);
df1_da2 = @(a1, a2) l3*cos(a1 + a2 - 2*pi/3);

dU_da1 = @(a1) k*(a1 - rest_alpha1);
dU_da2 = @(a2) k*(a2 - rest_alpha2);

% Function to solve for critical points
crit_system = @(x) [
    f1(x(1), x(2)); 
    dU_da1(x(1)) - x(3)*df1_da1(x(1), x(2));
    dU_da2(x(2)) - x(3)*df1_da2(x(1), x(2))
];

% Initial guesses for [α₁, α₂, λ₁]
init_guesses = [
    rest_alpha1, rest_alpha2, 0;       % Rest position
    rest_alpha1 + 0.5, rest_alpha2 - 0.5, 0; % Perturbed
    rest_alpha1 - 0.5, rest_alpha2 + 0.5, 0  % Perturbed
];

extrema_delta = [];
extrema_U = [];
extrema_type = {}; % 'min' or 'max'

tol = 1e-3; % Tolerance for duplicate detection

for i = 1:size(init_guesses, 1)
    [sol, ~, exitflag] = fsolve(crit_system, init_guesses(i,:)', options);
    
    if exitflag > 0
        a1_sol = sol(1);
        a2_sol = sol(2);
        delta_sol = -f2(a1_sol, a2_sol); % δ = -f₂
        
        % Check δ > 0 and not duplicate
        if delta_sol > 0 && ~any(abs(extrema_delta - delta_sol) < tol)
            U_sol = U(a1_sol, a2_sol);
            
            % Classify as min/max using curvature
            h = 1e-3;
            U_plus = interp1(delta_vec, U_vec, delta_sol + h, 'spline', 'extrap');
            U_minus = interp1(delta_vec, U_vec, delta_sol - h, 'spline', 'extrap');
            curvature = (U_plus - 2*U_sol + U_minus) / h^2;
            
            extrema_delta(end+1) = delta_sol;
            extrema_U(end+1) = U_sol;
            extrema_type{end+1} = 'Local Maximum';
            if curvature > 0
                extrema_type{end} = 'Local Minimum';
            end
        end
    end
end

%% Plot Results
figure;
plot(delta_vec, U_vec, 'b-', 'LineWidth', 2);
hold on;
scatter(extrema_delta, extrema_U, 100, 'r', 'filled');
grid on;
xlabel('Displacement \delta');
ylabel('Energy U');
title('Energy-Displacement Curve (\delta > 0)');
legend('U(\delta)', 'Critical Points', 'Location', 'best');

% Annotate extrema types
for i = 1:length(extrema_delta)
    text(extrema_delta(i), extrema_U(i), extrema_type{i}, ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

%% Display Extrema Values
disp('===== Local Extrema =====');
for i = 1:length(extrema_delta)
    fprintf('%s at δ = %.4f, U = %.4f\n', ...
            extrema_type{i}, extrema_delta(i), extrema_U(i));
end
