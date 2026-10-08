function cfg = unit_config
%UNIT_CONFIG Unit geometry and target edge stretches for main_unit_energy.
cfg.edge_length=15;
cfg.flank_length=12.75;
cfg.flank_width=0.75;
cfg.ligament_width=0.225;
cfg.beta=[0 pi/40];
cfg.edge_stretches=[1.63 1.63 1.63]; % [lambda12 lambda23 lambda31]; change for anisotropy
cfg.make_plot=true;
end
