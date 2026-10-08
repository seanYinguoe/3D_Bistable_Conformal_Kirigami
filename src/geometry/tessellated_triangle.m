function [triangle_tessellation,diagnostics] = tessellated_triangle(f_out, i_out, params, v, beta, varargin)
% Find the optimised bistable unit for each triangular unit, and move it to
% cooresponding grid

%% Generate optimised bistable unit for each triangular grid
% Create standards bistable unit(upwards and downwards)
edgeLen = params(1);
l1 = params(2);
l4 = params(3);
t = params(4);

if nargin >= 6 && ~isempty(varargin{1})
    t = varargin{1};
end

if isscalar(t)
    t = repmat(t, size(f_out,1), 1);
else
    t = t(:);
    if numel(t) ~= size(f_out,1)
        error('t must be a scalar or an Nx1 vector with one value per unit.');
    end
end

%% Generate optimised bistable unit for each grid and move bistable unit to responding grid
triangle_tessellation = cell(size(f_out,1),1);
diagnostics=struct('exitflag',zeros(size(f_out,1),1),'max_constraint',zeros(size(f_out,1),1));
for i = 1:size(f_out,1)
    % Calculate the edge1, edge2, and edge3 of each triangular grid
    q1 = v(f_out(i,1),:);
    q2 = v(f_out(i,2),:);
    q3 = v(f_out(i,3),:);
    if i_out(i) == 0 % upwards triangle
        [triangle,~,~,d] = deform_triangle(q3,q1,q2,edgeLen,l1,l4,beta(i),t(i),i_out(i));
    else % downwards triangle
        [triangle,~,~,d] = deform_triangle(q1,q3,q2,edgeLen,l1,l4,beta(i),t(i),i_out(i));
    end
    triangle_tessellation{i} = triangle;
    diagnostics.exitflag(i)=d.exitflag;diagnostics.max_constraint(i)=d.max_constraint;
end
if any(diagnostics.exitflag<=0|diagnostics.max_constraint>1e-6)
    warning('kirigami:DeployedGeometry','Some deployed-unit solves need review; inspect geometry_diagnostics.');
end
end
