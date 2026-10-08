function result = demo_unit_energy(makePlot)
%DEMO_UNIT_ENERGY Illustrative isotropic deployment using the existing HBM.
% Stores raw energies; normalisation in the plot is only for comparison.
% Example parameters are not a calibration of the published experiments.
if nargin < 1, makePlot = true; end
root = setup_project;
assert(exist('fmincon','file')==2 && exist('fsolve','file')==2,'Optimization Toolbox is required.');
cfg = struct('edge_length',15,'flank_length',12.75,'flank_width',0.75, ...
    'ligament_width',0.15,'segments',8,'beta',[0 pi/40],'points',81);
strain = linspace(0,0.65,cfg.points);
energy = zeros(numel(strain),numel(cfg.beta));
diagnostics = cell(size(cfg.beta));
for k=1:numel(cfg.beta)
    [~,E,diagnostics{k}] = deform_triangle_isotropic(strain*cfg.edge_length, ...
        cfg.edge_length,cfg.flank_length,cfg.flank_width,cfg.beta(k),cfg.ligament_width,cfg.segments,false);
    energy(:,k) = E(:);
    assert(all(diagnostics{k}.exitflag>0), 'An energy solve did not converge.');
    assert(all(diagnostics{k}.max_constraint<1e-6), 'Energy-solve constraints are not satisfied.');
end
assert(all(isfinite(energy(:))) && all(energy(:)>=0),'Invalid HBM energy output.');
out = tempname(fullfile(root,'results')); mkdir(out);
result = struct('config',cfg,'strain',strain,'energy',energy, ...
    'diagnostics',{diagnostics},'output_dir',out,'matlab_version',version);
save(fullfile(out,'unit_energy.mat'),'result');
writematrix([strain(:),energy],fullfile(out,'unit_energy.csv'));
if makePlot
    fig=figure('Color','w','Position',[100 100 760 450]);
    plot(strain,energy./max(energy,[],1),'LineWidth',1.8);
    xlabel('Isotropic edge strain'); ylabel('Energy / maximum of each curve');
    legend('beta = 0','beta = pi/40','Location','northwest');
    title('Illustrative HBM deployment sweep'); box off;
    exportgraphics(fig,fullfile(out,'unit_energy.png'),'Resolution',180);
end
fprintf('Unit energy example saved to %s\n',out);
end
