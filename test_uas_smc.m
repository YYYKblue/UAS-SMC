function tests = test_uas_smc
%TEST_UAS_SMC Regression and acceptance tests for the UAS-SMC simulation.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(mfilename('fullpath'));
addpath(root);
testCase.TestData.root = root;
end

function testRequiredImplementationFilesExist(testCase)
required = { ...
    'uas_smc_controller.m', ...
    'run_uas_smc_sim.m', ...
    'optimize_uas_smc.m'};
for idx = 1:numel(required)
    verifyEqual(testCase, exist(fullfile(testCase.TestData.root, required{idx}), 'file'), 2, ...
        sprintf('Missing required implementation file: %s', required{idx}));
end
end

function testControllerOriginIsEquilibrium(testCase)
[model, ctrl] = controllerFixture();
[u, diagnostic] = uas_smc_controller(zeros(4, 1), 0, model, ctrl);

verifyEqual(testCase, u, 0, 'AbsTol', 1e-14);
verifyEqual(testCase, diagnostic.phi, 0, 'AbsTol', 1e-14);
verifyEqual(testCase, diagnostic.N, 0, 'AbsTol', 1e-14);
verifyEqual(testCase, diagnostic.sigma, 0, 'AbsTol', 1e-14);
end

function testControllerIsOddAndReachingLawIsNonnegative(testCase)
[model, ctrl] = controllerFixture();
x = [0.01; -0.02; deg2rad(1.5); 0.03];
eta = -0.004;
[uPositive, diagnosticPositive] = uas_smc_controller(x, eta, model, ctrl);
[uNegative, diagnosticNegative] = uas_smc_controller(-x, -eta, model, ctrl);

verifyEqual(testCase, uNegative, -uPositive, 'AbsTol', 1e-10);
verifyGreaterThanOrEqual(testCase, diagnosticPositive.N, -1e-10);
verifyGreaterThanOrEqual(testCase, diagnosticNegative.N, -1e-10);

sigmaValues = linspace(-0.1, 0.1, 21);
etaValues = linspace(-0.05, 0.05, 17);
for sigma = sigmaValues
    xForSigma = model.S' * (sigma / (model.S * model.S'));
    for etaGrid = etaValues
        [~, diagnostic] = uas_smc_controller(xForSigma, etaGrid, model, ctrl);
        verifyGreaterThanOrEqual(testCase, diagnostic.N, -1e-9);
        verifyTrue(testCase, isfinite(diagnostic.phi));
    end
end
end

function testSwitchingLawIsContinuousAtBoundaries(testCase)
[model, ctrl] = controllerFixture();
delta = 1e-8;
eta = 0.01;
boundaries = [-ctrl.xi1 * eta, -ctrl.xi2 * eta];

for boundary = boundaries
    xMinus = model.S' * ((boundary - delta) / (model.S * model.S'));
    xPlus = model.S' * ((boundary + delta) / (model.S * model.S'));
    [~, diagnosticMinus] = uas_smc_controller(xMinus, eta, model, ctrl);
    [~, diagnosticPlus] = uas_smc_controller(xPlus, eta, model, ctrl);
    verifyLessThan(testCase, abs(diagnosticPlus.phi - diagnosticMinus.phi), 1e-6);
end
end

function testBaselineModelAndNominalSimulation(testCase)
simulation_options = struct( ...
    'useOptimized', false, ...
    'makePlots', false, ...
    'saveResults', false, ...
    'verbose', false, ...
    'Tend', 8, ...
    'x0', [0; 0; deg2rad(5); 0]); %#ok<NASGU>
run(fullfile(testCase.TestData.root, 'run_uas_smc_sim.m'));

expectedA = [0, 1, 0, 0; 0, -135.3752, -1.9188, 0; ...
             0, 0, 0, 1; 0, 380.1344, 32.9066, 0];
expectedB = [0; 10.5032; 0; -29.4932];
verifyEqual(testCase, sim_results.model.A, expectedA, 'AbsTol', 1e-4);
verifyEqual(testCase, sim_results.model.B, expectedB, 'AbsTol', 1e-4);
verifyEqual(testCase, rank(sim_results.model.Co), 4);
verifyGreaterThan(testCase, max(real(eig(sim_results.model.A))), 0);
verifyEqual(testCase, sim_results.selected.design.SB, 1, 'AbsTol', 1e-10);
verifyEqual(testCase, sort(sim_results.selected.design.internalPoles), [-4 -3 -2], 'AbsTol', 1e-8);

verifyEqual(testCase, numel(sim_results.selected.t), 1601);
verifyEqual(testCase, diff(sim_results.selected.t), ...
    0.005 * ones(1600, 1), 'AbsTol', 1e-12);
verifyLessThan(testCase, max(abs(sim_results.selected.state(end-199:end, 1))), 1e-3);
verifyLessThan(testCase, max(abs(sim_results.selected.state(end-199:end, 3))), deg2rad(0.1));
verifyLessThanOrEqual(testCase, sim_results.selected.metrics.peakVoltage, 12);
verifyLessThanOrEqual(testCase, sim_results.selected.metrics.peakAngleDeg, 10);
verifyGreaterThanOrEqual(testCase, min(sim_results.selected.diagnostic.N), -1e-9);
end

function testSymmetricAndLinearRangeScenarios(testCase)
scenarios = [5, -5, 8, -8];
terminalStates = zeros(4, 4);
for idx = 1:numel(scenarios)
    simulation_options = struct( ...
        'useOptimized', false, ...
        'makePlots', false, ...
        'saveResults', false, ...
        'verbose', false, ...
        'Tend', 8, ...
        'x0', [0; 0; deg2rad(scenarios(idx)); 0]); %#ok<NASGU>
    run(fullfile(testCase.TestData.root, 'run_uas_smc_sim.m'));
    response = sim_results.selected;
    terminalStates(idx, :) = response.state(end, :);
    verifyLessThan(testCase, max(abs(response.state(end-199:end, 1))), 1e-3);
    verifyLessThan(testCase, max(abs(response.state(end-199:end, 3))), deg2rad(0.1));
    verifyLessThanOrEqual(testCase, response.metrics.peakPosition, 0.20);
    verifyLessThanOrEqual(testCase, response.metrics.peakVoltage, 12);
    verifyLessThanOrEqual(testCase, response.metrics.peakAngleDeg, 10);
end
verifyEqual(testCase, terminalStates(2, :), -terminalStates(1, :), 'AbsTol', 1e-8);
verifyEqual(testCase, terminalStates(4, :), -terminalStates(3, :), 'AbsTol', 1e-8);
end

function testQuickOptimizerPreservesConstraintsAndBaseline(testCase)
outputDir = tempname;
optimizer_options = struct( ...
    'particles', 5, ...
    'iterations', 2, ...
    'trainingRandomCount', 2, ...
    'validationCount', 4, ...
    'localMaxRounds', 2, ...
    'saveResults', false, ...
    'verbose', false, ...
    'outputDir', outputDir); %#ok<NASGU>
run(fullfile(testCase.TestData.root, 'optimize_uas_smc.m'));

q = optimization_report.best.q;
verifySize(testCase, q, [1, 8]);
verifyGreaterThan(testCase, q(1), 0);
verifyGreaterThan(testCase, q(2), 0);
verifyGreaterThan(testCase, q(3), 0);
verifyGreaterThan(testCase, q(4), 0);
verifyGreaterThan(testCase, q(5), 0);
verifyGreaterThan(testCase, q(6), 0);
verifyGreaterThan(testCase, q(7), 0);
verifyGreaterThan(testCase, q(8), 0);
verifyGreaterThan(testCase, optimization_report.best.ctrl.xi1, optimization_report.best.ctrl.xi2);
verifyLessThan(testCase, optimization_report.best.ctrl.a, 0);
verifyGreaterThan(testCase, optimization_report.best.ctrl.b, 0);
verifyLessThanOrEqual(testCase, optimization_report.best.score, ...
    optimization_report.baseline.score + 1e-12);
verifyEqual(testCase, optimization_report.metadata.seed, 20260717);
verifyTrue(testCase, all(isfinite(optimization_report.history.bestScore)));
end

function testPhasePlaneAndReachingLawArtifacts(testCase)
outputDir = tempname;
simulation_options = struct( ...
    'useOptimized', false, ...
    'makePlots', false, ...
    'saveResults', true, ...
    'verbose', false, ...
    'Tend', 1, ...
    'outputDir', outputDir); %#ok<NASGU>
run(fullfile(testCase.TestData.root, 'run_uas_smc_sim.m'));

response = sim_results.selected;
verifyTrue(testCase, isfield(response.diagnostic, 'N_dot'));
verifySize(testCase, response.diagnostic.N_dot, size(response.diagnostic.N));
verifyTrue(testCase, all(isfinite(response.diagnostic.N_dot)));
verifyEqual(testCase, response.diagnostic.N_dot, ...
    gradient(response.diagnostic.N, sim_results.options.Ts), 'AbsTol', 1e-12);

required = {'uas_smc_phase_plane.png', 'uas_smc_reaching_law.png'};
for idx = 1:numel(required)
    info = dir(fullfile(outputDir, required{idx}));
    verifyNotEmpty(testCase, info, sprintf('Missing result artifact: %s', required{idx}));
    verifyGreaterThan(testCase, info.bytes, 0);
end
end

function testSavedFullOptimizationAndArtifacts(testCase)
resultsDir = fullfile(testCase.TestData.root, 'results');
required = { ...
    'uas_smc_best_params.mat', ...
    'optimization_history.mat', ...
    'optimization_convergence.png', ...
    'baseline_vs_optimized.png', ...
    'uas_smc_response.png', ...
    'uas_smc_diagnostics.png', ...
    'uas_smc_phase_plane.png', ...
    'uas_smc_reaching_law.png', ...
    'uas_smc_results.mat'};
for idx = 1:numel(required)
    info = dir(fullfile(resultsDir, required{idx}));
    verifyNotEmpty(testCase, info, sprintf('Missing result artifact: %s', required{idx}));
    verifyGreaterThan(testCase, info.bytes, 0);
end

optimized = load(fullfile(resultsDir, 'uas_smc_best_params.mat'));
report = optimized.optimization_report;
verifyEqual(testCase, report.metadata.particles, 40);
verifyEqual(testCase, report.metadata.iterations, 60);
verifyEqual(testCase, report.metadata.trainingScenarioCount, 16);
verifyEqual(testCase, report.metadata.validationScenarioCount, 100);
verifyLessThanOrEqual(testCase, report.best.score, report.baseline.score);
verifyTrue(testCase, all(report.best.evaluation.feasible));
verifyGreaterThanOrEqual(testCase, report.validation.passRate, 0.95);

validation = report.validation.evaluation;
verifyTrue(testCase, all(validation.allFinite));
verifyLessThanOrEqual(testCase, max(validation.peakAngleDeg), 10 + 1e-9);
verifyLessThanOrEqual(testCase, max(validation.peakPosition), 0.20 + 1e-9);
verifyLessThanOrEqual(testCase, max(validation.peakVoltage), 12 + 1e-9);
verifyGreaterThanOrEqual(testCase, min(validation.minimumN), -1e-9);

simulation = load(fullfile(resultsDir, 'uas_smc_results.mat'));
metrics = simulation.sim_results.selected.metrics;
verifyTrue(testCase, metrics.feasible);
verifyEqual(testCase, simulation.sim_results.options.Ts, 0.005, 'AbsTol', 1e-14);
verifyLessThan(testCase, metrics.settlingTime, 8);
verifyLessThanOrEqual(testCase, metrics.peakVoltage, 12);
end

function [model, ctrl] = controllerFixture()
[A, B] = nominalMatrices();
Co = [B, A * B, A^2 * B, A^3 * B];
h = [0 0 0 1] / Co;
transform = [h; h * A; h * A^2; h * A^3];
S = [24 26 9 1] * transform;
S = S / (S * B);
model = struct('A', A, 'B', B, 'S', S, 'SB', S * B);
ctrl = struct( ...
    'xi1', 2, ...
    'xi2', 1, ...
    'a', -1, ...
    'b', 1, ...
    'k', 4, ...
    'alpha', 1, ...
    'beta', 1, ...
    'm_aux', 0.01, ...
    'u_max', 12, ...
    'tol', 1e-12);
end

function [A, B] = nominalMatrices()
m = 0.0923;
M = 0.1945;
l = 0.1950;
J = 0.0029;
g = 9.8;
R = 3.75;
r = 0.018;
Km = 0.183;
Ke = 0.232;
I = 7.083e-6;
Qeq = m * J + (J + m * l^2) * (M + I / r^2);
A22 = -(Km * Ke * (J + m * l^2)) / (Qeq * R * r^2);
A23 = -(m^2 * l^2 * g) / Qeq;
A42 = (m * l * Km * Ke) / (Qeq * R * r^2);
A43 = (m * g * l * (M + m + I / r^2)) / Qeq;
B21 = (Km * (J + m * l^2)) / (Qeq * R * r);
B41 = -(m * l * Km) / (Qeq * R * r);
A = [0 1 0 0; 0 A22 A23 0; 0 0 0 1; 0 A42 A43 0];
B = [0; B21; 0; B41];
end
