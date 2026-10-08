function tests = test_geometry
tests = functiontests(localfunctions);
end
function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
setup_project;
end
function testAreaScalingIsConsistent(testCase)
v=[0 0;1 0;0 1];f=[1 2 3];
testCase.verifyEqual(triangle_area_2D(v,f),0.5,'AbsTol',1e-14);
testCase.verifyEqual(triangle_area_3D([v zeros(3,1)]*2,f),2,'AbsTol',1e-14);
end
function testEdgeScalingMatchesUniformExpansion(testCase)
v=[0 0 0;1 0 0;0.5 sqrt(3)/2 0];f=[1 2 3];
s=calculate_scale_facs(v,1.5*v,f);
testCase.verifyEqual(s,1.5*ones(size(s)),'AbsTol',1e-12);
end
function testMeshAndUvConnectivityAgree(testCase)
root=setup_project;
for model={'quarter_dome','double_dome','hemisphere'}
    name=model{1};obj=readObj(fullfile(root,'data','meshes',name),[name '_flat.obj']);
    uv=vertice_sort(obj.vt,obj.f.v,obj.f.vt);
    testCase.verifyEqual(size(uv,1),size(obj.v,1));
    testCase.verifyTrue(all(isfinite(uv(:))));
    testCase.verifyGreaterThan(triangle_area_2D(uv,obj.f.v),0);
    testCase.verifyGreaterThan(triangle_area_3D(obj.v,obj.f.v),0);
end
end
function testLibraryHasRequiredColumns(testCase)
root=setup_project;
s=load(fullfile(root,'data','unit-library','anisotropy_filter_refine.mat'));
T=s.anisotropy_filter_refine;
testCase.verifyTrue(all(ismember({'a1','a2','a3','eps_bist','eta_val','beta'},T.Properties.VariableNames)));
testCase.verifyGreaterThan(height(T),0);
testCase.verifyEqual(T.a1+T.a2+T.a3,pi*ones(height(T),1),'AbsTol',1e-10);
end
function testTargetTrianglePreservesEdgeOrdering(testCase)
lambda = [1.2 1.4 1.3];
[q1, q2, q3, valid] = scale_facs_to_q(lambda, 15);
testCase.verifyTrue(valid);
actual = [norm(q1-q2), norm(q2-q3), norm(q3-q1)];
testCase.verifyEqual(actual, 15*lambda, 'AbsTol', 1e-12);
end
function testImpossibleTargetTriangleIsRejected(testCase)
[~, ~, ~, valid, reason] = scale_facs_to_q([1 1 3], 15);
testCase.verifyFalse(valid);
testCase.verifyEqual(reason, 'invalid_target_triangle');
end
function testRimAssignmentIsExplicit(testCase)
angles = repmat(pi/3, 6, 1);
T = table(angles, angles, angles, linspace(.3,.6,6)', .8*ones(6,1), ...
    linspace(0,pi/20,6)', 'VariableNames', {'a1','a2','a3','eps_bist','eta_val','beta'});
query = [pi/3 pi/3 pi/3 .45; pi/3 pi/3 pi/3 2];
[beta, eta, info] = assign_opt_beta(query, T, struct('rimMask', [true;false]));
testCase.verifyEqual(beta(1), pi/12);
testCase.verifyTrue(isnan(eta(1)));
testCase.verifyEqual(info.flag(1), "rim_undeployed");
% A nearest library match must retain a large residual, not imply success.
testCase.verifyGreaterThan(info.beta_error(2), 1);
end
function testTargetRescalingPreservesStretchRatios(testCase)
v = [0 0 0; 1 0 0; .5 sqrt(3)/2 0];
f = [1 2 3];
target = [0 0 0; 1.2 0 .1; .4 1 .2];
before = calculate_scale_facs(v, target, f);
[~, after, info] = rescale_target_edges(v, f, struct(), target, 1.4, struct());
testCase.verifyEqual(after, 1.4*before, 'AbsTol', 1e-12);
testCase.verifyEqual(info.edge_scale_check, [1.4 1.4], 'AbsTol', 1e-12);
end
