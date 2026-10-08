function cfg = surface_config
%SURFACE_CONFIG Input and geometry settings for main_surface.
root=fileparts(fileparts(mfilename('fullpath')));
cfg.mesh=fullfile(root,'data','meshes','quarter_dome','quarter_dome_flat.obj');
cfg.library=fullfile(root,'data','unit-library','anisotropy_filter_refine.mat');
cfg.edge_length=16;
cfg.flank_length_ratio=0.85;
cfg.flank_width_ratio=0.05;
cfg.ligament_width_ratios=[0.015 0.025];
cfg.thickness_threshold=0.4;
cfg.library_eta_min=0.10;
cfg.assignment=struct('kNN',80,'nBetaEval',61,'strainTol',0.015);
cfg.build_deployed=true;
cfg.export_svg=true;
cfg.fillet_radius=0.30; % same coordinate units as the mesh, not an automatic kerf allowance
cfg.make_plot=true;
cfg.evaluate_energy=false; % optional independent-unit HBM sum; can be expensive
cfg.energy=struct('points',150,'segments',8,'use_parallel',false);
end
