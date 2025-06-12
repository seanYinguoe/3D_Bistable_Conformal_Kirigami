function plotgrid(f, v, varargin)
%   Visualize 2D triangulation with face numbers
%   PLOTGRID(F, V) plots mesh with face numbers using:
%       F - Mx3 connectivity matrix
%       V - Nx2 vertex coordinates
%
%   PLOTGRID(F, V, Name,Value) specifies plot parameters:
%       'EdgeColor'    - Edge color (default: 'b')
%       'LineWidth'    - Edge line width (default: 1)
%       'FontSize'     - Face number size (default: 9)
%       'TextColor'    - Text color (default: 'k')
%       'Background'   - Text background (default: [1 1 1 0.7])

% Input parsing
p = inputParser;
addParameter(p, 'EdgeColor', 'b', @(x) ischar(x) || isvector(x));
addParameter(p, 'LineWidth', 1, @isscalar);
addParameter(p, 'FontSize', 12, @isscalar);
addParameter(p, 'TextColor', 'k', @(x) ischar(x) || isvector(x));
addParameter(p, 'Background', [1 1 1 0.7], @(x) validateattributes(x, {'numeric'}, {'vector','numel',4}));
parse(p, varargin{:});

% Create patch plot
patch('Vertices', v, 'Faces', f,...
      'FaceColor', 'none',...
      'EdgeColor', p.Results.EdgeColor,...
      'LineWidth', p.Results.LineWidth);
hold on;

% Calculate centers using triangulation
TR = triangulation(f, v);
centers = incenter(TR);

% Add face numbers
text(centers(:,1), centers(:,2), string(1:size(f,1)),...
    'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'middle',...
    'FontSize', p.Results.FontSize,...
    'Color', p.Results.TextColor,...
    'BackgroundColor', p.Results.Background,...
    'Margin', 0.5);

hold off;
axis equal;
axis off
box on;
grid off;
end
