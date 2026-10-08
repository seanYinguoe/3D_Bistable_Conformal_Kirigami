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
