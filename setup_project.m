function root = setup_project
%SETUP_PROJECT Configure this repository for the current MATLAB session only.
% Run from the repository root; does not change MATLAB's saved path.
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'config'), genpath(fullfile(root, 'src')), ...
    fullfile(root, 'third_party', 'Objread'));
if ~isfolder(fullfile(root, 'results'))
    mkdir(fullfile(root, 'results'));
end
end
