function out = model
%
% test.m
%
% Model exported on May 6 2025, 12:29 by COMSOL 6.2.0.290.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath(pwd);

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').physics.create('shell', 'Shell', 'geom1');

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');
model.study('std1').feature('stat').setSolveFor('/physics/shell', true);

model.component('comp1').geom('geom1').run;

model.param.set('tk', '0.3[m]');

model.component('comp1').physics('shell').feature('to1').set('d', 'tk');

out = model;
