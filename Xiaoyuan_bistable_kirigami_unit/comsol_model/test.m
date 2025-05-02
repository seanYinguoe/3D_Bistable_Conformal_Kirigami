function out = model
%
% test.m
%
% Model exported on May 1 2025, 16:45 by COMSOL 6.2.0.290.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/Users/sean/Desktop/Project 2/3D_Bistable_Conformal_Kirigami/Xiaoyuan_bistable_kirigami_unit/comsol_model');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').geom.create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').set('table', [0 0; 0.5 0.5; 1 0]);
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').physics.create('shell', 'Shell', 'geom1');

model.component('comp1').view('view2').axis.set('xmin', -0.05000004172325134);
model.component('comp1').view('view2').axis.set('xmax', 1.0499999523162842);
model.component('comp1').view('view2').axis.set('ymin', -0.10835176706314087);
model.component('comp1').view('view2').axis.set('ymax', 0.6083517670631409);

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');

model.component('comp1').geom('geom1').feature('wp1').geom.feature.duplicate('pol2', 'pol1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol2').setIndex('table', 0.3, 1, 0);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol2').setIndex('table', 0.3, 1, 1);
model.component('comp1').geom('geom1').feature('wp1').geom.run('pol2');

out = model;
