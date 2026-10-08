function root = setup_project
%SETUP_PROJECT Add this paper's paths without changing the saved MATLAB path.
root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(genpath(fullfile(root,'src')));
addpath(fullfile(root,'third_party','Objread'));
addpath(fullfile(root,'examples'));
if ~isfolder(fullfile(root,'results')), mkdir(fullfile(root,'results')); end
end
