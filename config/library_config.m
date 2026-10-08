function cfg = library_config
%LIBRARY_CONFIG Dense research sweep; use a small grid when checking the setup.
cfg.alpha3=10*pi/36:pi/144:12*pi/36;
cfg.alpha2=10*pi/36:pi/144:13*pi/36;
cfg.beta=linspace(0,pi/20,20);
cfg.edge_length=15;
cfg.scale=1.63;
cfg.flank_length_ratio=0.80; % current sweep setting; historical supplied table used 0.85
cfg.flank_width_ratio=0.05;
cfg.ligament_width_ratio=0.015;
cfg.use_parallel=false;
cfg.postprocess=true;
cfg.refine=true;
end
