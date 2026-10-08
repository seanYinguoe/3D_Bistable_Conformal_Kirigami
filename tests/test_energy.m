function tests = test_energy
% Cross-check the three-ligament model against its isotropic reduction.
tests = functiontests(localfunctions);
end
function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
setup_project;
end
function testIsotropicLimitMatchesReducedModel(testCase)
edge = 15;
beta = pi/40;
strain = linspace(0, .55, 41);
[~, expected, diagnostics] = deform_triangle_isotropic( ...
    strain*edge, edge, 12.75, .75, beta, .225, 8, false);
[q1, q2, q3] = scale_facs_to_q(1.55*ones(1,3), edge);
actual = deform_triangle_anisotropic(q1, q2, q3, edge, 12.75, .75, beta, .225, 41, 8, false);
testCase.verifyGreaterThan(diagnostics.exitflag, 0);
testCase.verifyLessThan(diagnostics.max_constraint, 1e-6);
testCase.verifyLessThan(max(abs(actual-expected))/max(expected), 1e-5);
end
function testInvalidUnitCannotDisappearFromEnergySum(testCase)
options = struct('points', 5, 'segments', 8, 'use_parallel', false);
[total, units, info] = calculate_global_energy([15 12.75 .75 .225], pi/40, .225, ...
    [1.04 1.04 1.04; 1 1 3], false, options);
testCase.verifyTrue(all(isfinite(units(1,:))));
testCase.verifyTrue(all(isnan(total)));
testCase.verifyEqual(info.N_invalid_units, 1);
testCase.verifyEqual(info.invalid_unit_index, 2);
end
