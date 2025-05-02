% Connect matlab and Comsol
% Input: geometry model of optimised kirigami tessellation
clc;
clear;
import com.comsol.model.*
import com.comsol.model.util.*
%% Geometry modelling
% Load data
load('tessellation.mat');
% Define connecivity of filaments, voids, flanks, Innertriangle
f_void = [1  2  3  4  5  6;
 7  8  9 10 11 12;
13 14 15 16 17 18];
f_flank = [19 20 21 22;
23 24 25 26;
27 28 29 30];
f_filament = [31 32 33 34;
35 36 37 38;
39 40 41 42];
f_inner = [43 44 45];

model = ModelUtil.create('Model');  % Define a model
comp1 = model.component.create('comp1', true);
geom1 = comp1.geom.create('geom1', 3);  % Correct geometry reference
wp1 = geom1.feature.create('wp1', 'WorkPlane');  % Use 'geom1' not 'geom1'
wp1.set('planetype', 'quick');
wp1.set('quickplane', 'xy');
% Build different parts of tessellation
for i = 1:size(tessellation,1)
    for j = 1:size(f_flank,1)
    % Build flanks
    name = strcat('gflanks', num2str(i), num2str(j));
    polygon = wp1.geom.feature.create(name, 'Polygon');
    x = tessellation{i}(f_flank(j,:), 1);
    y = tessellation{i}(f_flank(j,:), 2);
    polygon.set('x', x);
    polygon.set('y', y);
    end
    for j = 1:size(f_filament,1)
    % Build filament
    name = strcat('gfilament', num2str(i), num2str(j));
    polygon = wp1.geom.feature.create(name, 'Polygon');
    x = tessellation{i}(f_filament(j,:),1);
    y = tessellation{i}(f_filament(j,:),2);
    polygon.set('x', x);
    polygon.set('y', y);
    end
    for j = 1:size(f_inner,1)
    % Build inner triangle
    name = strcat('ginner', num2str(i), num2str(j));
    polygon = wp1.geom.feature.create(name, 'Polygon');
    x = tessellation{i}(f_inner(j,:),1);
    y = tessellation{i}(f_inner(j,:),2);
    polygon.set('x', x);
    polygon.set('y', y);
    end
end
model.component('comp1').geom('geom1').runAll; % build model
mphgeom(model, 'geom1'); % show the geometry results

%% Boundary condition

%% Mesh

%% Study

%% Post processing
