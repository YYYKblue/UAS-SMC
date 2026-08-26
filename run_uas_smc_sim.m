%RUN_UAS_SMC_SIM Simulate the 200 Hz linear inverted-pendulum UAS-SMC.
% Optional caller variable: simulation_options (structure).  The script
% leaves sim_results in the caller workspace and saves plots/data by default.

scriptRoot = fileparts(mfilename('fullpath'));
defaultOptions = struct( ...
    'Ts', 0.005, ...
    'Tend', 8, ...
    'x0', [0; 0; deg2rad(5); 0], ...
    'eta0', 0, ...
    'useOptimized', true, ...
    'ctrl', [], ...
    'modelParams', struct(), ...
    'makePlots', true, ...
    'saveResults', true, ...
    'verbose', true, ...
    'figureVisible', 'off', ...
    'outputDir', fullfile(scriptRoot, 'results'));

if ~exist('simulation_options', 'var') || isempty(simulation_options)
    simulation_options = defaultOptions;
else
    simulation_options = localMergeStruct(defaultOptions, simulation_options);
end

params = localMergeStruct(localNominalParameters(), simulation_options.modelParams);
model = localBuildModel(params);
ctrlBaseline = localBaselineController();

if ~isempty(simulation_options.ctrl)
    ctrlSelected = localNormalizeController(simulation_options.ctrl);
    selectedLabel = 'caller-supplied';
elseif simulation_options.useOptimized
    parameterFile = fullfile(simulation_options.outputDir, 'uas_smc_best_params.mat');
    if exist(parameterFile, 'file') ~= 2
        error('uas_smc:MissingOptimizedParameters', ...
            ['Optimized parameters are missing. Run optimize_uas_smc.m first, ' ...
             'or set simulation_options.useOptimized=false.']);
    end
    loaded = load(parameterFile, 'best_ctrl');
    ctrlSelected = localNormalizeController(loaded.best_ctrl);
    selectedLabel = 'optimized';
else
    ctrlSelected = ctrlBaseline;
    selectedLabel = 'baseline';
end

baselineResponse = localSimulate(model, ctrlBaseline, simulation_options);
if isequaln(ctrlSelected, ctrlBaseline)
    selectedResponse = baselineResponse;
else
    selectedResponse = localSimulate(model, ctrlSelected, simulation_options);
end

sim_results = struct();
sim_results.model = model;
sim_results.options = simulation_options;
sim_results.baseline = baselineResponse;
sim_results.selected = selectedResponse;
sim_results.selectedLabel = selectedLabel;
sim_results.ctrlBaseline = ctrlBaseline;
sim_results.ctrlSelected = ctrlSelected;

if simulation_options.makePlots || simulation_options.saveResults
    if exist(simulation_options.outputDir, 'dir') ~= 7
        mkdir(simulation_options.outputDir);
    end
    figures = localCreatePlots(baselineResponse, selectedResponse, selectedLabel, ...
        simulation_options.figureVisible);
    if simulation_options.saveResults
        exportgraphics(figures.response, ...
            fullfile(simulation_options.outputDir, 'uas_smc_response.png'), 'Resolution', 180);
        exportgraphics(figures.diagnostics, ...
            fullfile(simulation_options.outputDir, 'uas_smc_diagnostics.png'), 'Resolution', 180);
        exportgraphics(figures.comparison, ...
            fullfile(simulation_options.outputDir, 'baseline_vs_optimized.png'), 'Resolution', 180);
        exportgraphics(figures.phasePlane, ...
            fullfile(simulation_options.outputDir, 'uas_smc_phase_plane.png'), 'Resolution', 180);
        exportgraphics(figures.reachingLaw, ...
            fullfile(simulation_options.outputDir, 'uas_smc_reaching_law.png'), 'Resolution', 180);
    end
    if ~simulation_options.makePlots || strcmpi(simulation_options.figureVisible, 'off')
        close(figures.response);
        close(figures.diagnostics);
        close(figures.comparison);
        close(figures.phasePlane);
        close(figures.reachingLaw);
    end
end

if simulation_options.saveResults
    if exist(simulation_options.outputDir, 'dir') ~= 7
        mkdir(simulation_options.outputDir);
    end
    save(fullfile(simulation_options.outputDir, 'uas_smc_results.mat'), 'sim_results');
end

if simulation_options.verbose
    metrics = selectedResponse.metrics;
    fprintf('\nUAS-SMC simulation (%s)\n', selectedLabel);
    fprintf('  Ts / frequency       : %.6f s / %.1f Hz\n', ...
        simulation_options.Ts, 1 / simulation_options.Ts);
    fprintf('  Internal poles       : [%s]\n', num2str(selectedResponse.design.internalPoles, ' %.4g'));
    fprintf('  xi1, xi2, a, b, k    : %.5g, %.5g, %.5g, %.5g, %.5g\n', ...
        ctrlSelected.xi1, ctrlSelected.xi2, ctrlSelected.a, ctrlSelected.b, ctrlSelected.k);
    fprintf('  Final x / theta      : %.6g m / %.6g deg\n', ...
        metrics.finalPosition, metrics.finalAngleDeg);
    fprintf('  Peak |x| / |theta|   : %.6g m / %.6g deg\n', ...
        metrics.peakPosition, metrics.peakAngleDeg);
    fprintf('  Peak voltage         : %.6g V\n', metrics.peakVoltage);
    fprintf('  Settling time        : %.6g s\n', metrics.settlingTime);
    fprintf('  Feasible             : %d\n\n', metrics.feasible);
end

function params = localNominalParameters()
params = struct( ...
    'm', 0.0923, ...
    'M', 0.1945, ...
    'l', 0.1950, ...
    'J', 0.0029, ...
    'g', 9.8, ...
    'R', 3.75, ...
    'r', 0.018, ...
    'Km', 0.183, ...
    'Ke', 0.232, ...
    'I', 7.083e-6);
end

function model = localBuildModel(p)
Qeq = p.m * p.J + (p.J + p.m * p.l^2) * (p.M + p.I / p.r^2);
A22 = -(p.Km * p.Ke * (p.J + p.m * p.l^2)) / (Qeq * p.R * p.r^2);
A23 = -(p.m^2 * p.l^2 * p.g) / Qeq;
A42 = (p.m * p.l * p.Km * p.Ke) / (Qeq * p.R * p.r^2);
A43 = (p.m * p.g * p.l * (p.M + p.m + p.I / p.r^2)) / Qeq;
B21 = (p.Km * (p.J + p.m * p.l^2)) / (Qeq * p.R * p.r);
B41 = -(p.m * p.l * p.Km) / (Qeq * p.R * p.r);
A = [0 1 0 0; 0 A22 A23 0; 0 0 0 1; 0 A42 A43 0];
B = [0; B21; 0; B41];
Co = [B, A * B, A^2 * B, A^3 * B];
if rank(Co) ~= 4
    error('uas_smc:UncontrollableModel', 'The supplied inverted-pendulum model is not controllable.');
end
model = struct('params', p, 'Qeq', Qeq, 'A', A, 'B', B, 'Co', Co);
end

function ctrl = localBaselineController()
ctrl = struct( ...
    'poles', [2 3 4], ...
    'xi1', 2, ...
    'xi2', 1, ...
    'alpha', 1, ...
    'beta', 1, ...
    'a', -1, ...
    'b', 1, ...
    'k', 4, ...
    'm_aux', 0.01, ...
    'u_max', 12, ...
    'tol', 1e-12);
end

function ctrl = localNormalizeController(ctrl)
required = {'poles', 'xi1', 'xi2', 'alpha', 'beta', 'k'};
for idx = 1:numel(required)
    if ~isfield(ctrl, required{idx})
        error('uas_smc:IncompleteController', 'Missing controller field: %s', required{idx});
    end
end
ctrl.a = -ctrl.alpha;
ctrl.b = ctrl.beta;
if ~isfield(ctrl, 'm_aux'), ctrl.m_aux = 0.01; end
if ~isfield(ctrl, 'u_max'), ctrl.u_max = 12; end
if ~isfield(ctrl, 'tol'), ctrl.tol = 1e-12; end
end

function design = localBuildDesign(model, ctrl)
poles = sort(ctrl.poles(:).');
if numel(poles) ~= 3 || any(poles <= 0) || any(diff(poles) <= 0)
    error('uas_smc:InvalidInternalPoles', 'Require three distinct, increasing positive pole magnitudes.');
end
h = [0 0 0 1] / model.Co;
transform = [h; h * model.A; h * model.A^2; h * model.A^3];
polynomial = poly(-poles);
surfaceCoefficients = [polynomial(4), polynomial(3), polynomial(2), 1];
S = surfaceCoefficients * transform;
S = S / (S * model.B);
design = struct( ...
    'S', S, ...
    'SB', S * model.B, ...
    'transform', transform, ...
    'surfaceCoefficients', surfaceCoefficients, ...
    'internalPoles', -poles);
end

function response = localSimulate(model, ctrl, options)
ctrl = localNormalizeController(ctrl);
design = localBuildDesign(model, ctrl);
controllerModel = struct('A', model.A, 'B', model.B, 'S', design.S, 'SB', design.SB);

Aaug = [model.A, zeros(4, 1); design.S, 0];
Baug = [model.B; 0];
transition = expm([Aaug, Baug; zeros(1, 6)] * options.Ts);
Ad = transition(1:5, 1:5);
Bd = transition(1:5, 6);

t = (0:options.Ts:options.Tend).';
sampleCount = numel(t);
z = zeros(5, sampleCount);
z(:, 1) = [options.x0(:); options.eta0];
u = zeros(sampleCount, 1);
diagnostic = struct();
fields = {'sigma', 's1', 's2', 'epsilon', 'phi', 'N', 'h', ...
    'region', 'omega1', 'omega2', 'u_unsat', 'saturated'};
for fieldIndex = 1:numel(fields)
    diagnostic.(fields{fieldIndex}) = zeros(sampleCount, 1);
end

for sample = 1:sampleCount
    [u(sample), current] = uas_smc_controller(z(1:4, sample), z(5, sample), ...
        controllerModel, ctrl);
    for fieldIndex = 1:numel(fields)
        diagnostic.(fields{fieldIndex})(sample) = current.(fields{fieldIndex});
    end
    if sample < sampleCount
        z(:, sample + 1) = Ad * z(:, sample) + Bd * u(sample);
    end
end
diagnostic.N_dot = gradient(diagnostic.N, options.Ts);

response = struct();
response.t = t;
response.state = z(1:4, :).';
response.eta = z(5, :).';
response.u = u;
response.diagnostic = diagnostic;
response.ctrl = ctrl;
response.design = design;
response.metrics = localMetrics(t, response.state, u, diagnostic, ctrl);
end

function metrics = localMetrics(t, state, u, diagnostic, ctrl)
position = state(:, 1);
angle = state(:, 3);
within = abs(position) < 1e-3 & abs(angle) < deg2rad(0.1);
stableSuffix = false(size(within));
stableSuffix(end) = within(end);
for idx = numel(within)-1:-1:1
    stableSuffix(idx) = within(idx) && stableSuffix(idx + 1);
end
settlingIndex = find(stableSuffix, 1, 'first');
if isempty(settlingIndex)
    settlingTime = NaN;
else
    settlingTime = t(settlingIndex);
end
tailCount = min(200, numel(t));
tail = (numel(t) - tailCount + 1):numel(t);
metrics = struct( ...
    'finalPosition', position(end), ...
    'finalAngleDeg', rad2deg(angle(end)), ...
    'peakPosition', max(abs(position)), ...
    'peakAngleDeg', max(abs(rad2deg(angle))), ...
    'peakVoltage', max(abs(u)), ...
    'peakRawVoltage', max(abs(diagnostic.u_unsat)), ...
    'rmsVoltage', sqrt(mean(u.^2)), ...
    'rmsDeltaVoltage', sqrt(mean(diff(u).^2)), ...
    'iaePosition', trapz(t, abs(position)), ...
    'iaeAngle', trapz(t, abs(angle)), ...
    'settlingTime', settlingTime, ...
    'tailPeakPosition', max(abs(position(tail))), ...
    'tailPeakAngleDeg', max(abs(rad2deg(angle(tail)))), ...
    'saturationFraction', mean(diagnostic.saturated), ...
    'minimumReachingLaw', min(diagnostic.N));
metrics.feasible = all(isfinite(state), 'all') ...
    && all(isfinite(u)) ...
    && metrics.peakPosition <= 0.20 ...
    && metrics.peakAngleDeg <= 10 ...
    && metrics.peakRawVoltage <= ctrl.u_max + 1e-9 ...
    && metrics.tailPeakPosition < 1e-3 ...
    && metrics.tailPeakAngleDeg < 0.1 ...
    && metrics.minimumReachingLaw >= -1e-9 ...
    && isfinite(metrics.settlingTime);
end

function figures = localCreatePlots(baseline, selected, selectedLabel, visibility)
figures.response = figure('Name', 'UAS-SMC response', 'Visible', visibility, ...
    'Color', 'w', 'Position', [100 100 1050 720]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; plot(selected.t, selected.state(:, 1), 'LineWidth', 1.4); grid on;
xlabel('Time (s)'); ylabel('x (m)'); title('Cart position');
nexttile; plot(selected.t, rad2deg(selected.state(:, 3)), 'LineWidth', 1.4); grid on;
xlabel('Time (s)'); ylabel('\theta (deg)'); title('Pendulum angle');
nexttile; plot(selected.t, selected.state(:, 2), 'LineWidth', 1.4); grid on;
xlabel('Time (s)'); ylabel('dx/dt (m/s)'); title('Cart velocity');
nexttile; plot(selected.t, selected.state(:, 4), 'LineWidth', 1.4); grid on;
xlabel('Time (s)'); ylabel('d\theta/dt (rad/s)'); title('Angular velocity');
nexttile([1 2]); plot(selected.t, selected.u, 'LineWidth', 1.4); grid on;
yline(12, '--r'); yline(-12, '--r'); xlabel('Time (s)'); ylabel('u_a (V)');
title(sprintf('Control voltage (%s)', selectedLabel));

figures.diagnostics = figure('Name', 'UAS-SMC diagnostics', 'Visible', visibility, ...
    'Color', 'w', 'Position', [120 120 1050 720]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; plot(selected.t, selected.diagnostic.sigma, 'LineWidth', 1.2); grid on;
xlabel('Time (s)'); ylabel('\sigma'); title('Virtual state');
nexttile; plot(selected.t, [selected.diagnostic.s1, selected.diagnostic.s2], 'LineWidth', 1.1); grid on;
xlabel('Time (s)'); ylabel('s'); legend('s_1', 's_2'); title('Switching surfaces');
nexttile; plot(selected.t, selected.diagnostic.h, 'LineWidth', 1.2); grid on;
xlabel('Time (s)'); ylabel('h'); title('Active auxiliary surface');
nexttile; plot(selected.t, selected.diagnostic.N, 'LineWidth', 1.2); grid on;
xlabel('Time (s)'); ylabel('N'); title('Reaching law');
nexttile; plot(selected.t, selected.diagnostic.epsilon, 'LineWidth', 1.2); grid on;
xlabel('Time (s)'); ylabel('\epsilon'); ylim([-0.05 1.05]); title('Continuous blending');
nexttile; stairs(selected.t, selected.diagnostic.region, 'LineWidth', 1.0); grid on;
xlabel('Time (s)'); ylabel('Region'); yticks(0:3); title('Active subspace');

figures.comparison = figure('Name', 'Baseline versus optimized', 'Visible', visibility, ...
    'Color', 'w', 'Position', [140 140 1050 650]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; plot(baseline.t, baseline.state(:, 1), '--', selected.t, selected.state(:, 1), '-', 'LineWidth', 1.25); grid on;
xlabel('Time (s)'); ylabel('x (m)'); legend('Baseline', selectedLabel); title('Cart position');
nexttile; plot(baseline.t, rad2deg(baseline.state(:, 3)), '--', selected.t, rad2deg(selected.state(:, 3)), '-', 'LineWidth', 1.25); grid on;
xlabel('Time (s)'); ylabel('\theta (deg)'); legend('Baseline', selectedLabel); title('Pendulum angle');
nexttile; plot(baseline.t, baseline.u, '--', selected.t, selected.u, '-', 'LineWidth', 1.25); grid on;
xlabel('Time (s)'); ylabel('u_a (V)'); legend('Baseline', selectedLabel); title('Voltage');
nexttile; plot(baseline.t, baseline.diagnostic.sigma, '--', selected.t, selected.diagnostic.sigma, '-', 'LineWidth', 1.25); grid on;
xlabel('Time (s)'); ylabel('\sigma'); legend('Baseline', selectedLabel); title('Virtual state');

figures.phasePlane = figure('Name', 'UAS-SMC phase plane', 'Visible', visibility, ...
    'Color', 'w', 'Position', [160 160 1000 720]);
ctrl = selected.ctrl;

% The positively invariant set is the intersection h_i >= 0 of the four
% auxiliary half-planes.  With the symmetric auxiliary surfaces this is
% equivalently |sigma + alpha*eta| <= m_aux and
% |sigma - beta*eta| <= m_aux.
auxiliaryCoordinates = ctrl.m_aux * [1 1; 1 -1; -1 -1; -1 1; 1 1];
invariantEta = (auxiliaryCoordinates(:, 1) - auxiliaryCoordinates(:, 2)) ...
    / (ctrl.alpha + ctrl.beta);
invariantSigma = (ctrl.beta * auxiliaryCoordinates(:, 1) ...
    + ctrl.alpha * auxiliaryCoordinates(:, 2)) / (ctrl.alpha + ctrl.beta);

etaMagnitude = max(abs(selected.eta));
etaExtent = max([1.15 * etaMagnitude, 1.15 * max(abs(invariantEta)), 1e-3]);
etaGrid = linspace(-etaExtent, etaExtent, 500);
hold on;
patch(invariantEta, invariantSigma, [0.72 0.88 0.68], ...
    'FaceAlpha', 0.55, 'EdgeColor', [0.15 0.50 0.15], ...
    'LineWidth', 1.5, 'DisplayName', 'Positive invariant set Q');
plot(etaGrid, -ctrl.xi1 * etaGrid, '-', 'Color', [0.85 0.15 0.15], ...
    'LineWidth', 1.8, 'DisplayName', 's_1 = 0');
plot(etaGrid, -ctrl.xi2 * etaGrid, '-', 'Color', [0.95 0.45 0.10], ...
    'LineWidth', 1.8, 'DisplayName', 's_2 = 0');
plot(etaGrid, -ctrl.alpha * etaGrid - ctrl.m_aux, ':', ...
    'Color', [0.10 0.45 0.80], 'LineWidth', 1.4, 'DisplayName', 'h_0 = 0');
plot(etaGrid, ctrl.beta * etaGrid + ctrl.m_aux, ':', ...
    'Color', [0.15 0.65 0.55], 'LineWidth', 1.4, 'DisplayName', 'h_1 = 0');
plot(etaGrid, ctrl.beta * etaGrid - ctrl.m_aux, ':', ...
    'Color', [0.40 0.30 0.75], 'LineWidth', 1.4, 'DisplayName', 'h_2 = 0');
plot(etaGrid, -ctrl.alpha * etaGrid + ctrl.m_aux, ':', ...
    'Color', [0.65 0.25 0.65], 'LineWidth', 1.4, 'DisplayName', 'h_3 = 0');
plot(selected.eta, selected.diagnostic.sigma, 'k-', 'LineWidth', 2.1, ...
    'DisplayName', 'System trajectory');
plot(selected.eta(1), selected.diagnostic.sigma(1), 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.20 0.70 0.20], ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Start');
plot(selected.eta(end), selected.diagnostic.sigma(end), 's', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.10 0.30 0.90], ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'End');

% Label the four subspaces using the same (s1, s2) sign convention as the
% controller: 0=(--), 1=(-+), 2=(+-), and 3=(++).
regionEtaMagnitude = 0.72 * etaExtent;
meanSwitchingSlope = 0.5 * (ctrl.xi1 + ctrl.xi2);
regionEta = [-regionEtaMagnitude; -regionEtaMagnitude; ...
    regionEtaMagnitude; regionEtaMagnitude];
regionSigma = [-0.25 * ctrl.xi2 * regionEtaMagnitude; ...
    meanSwitchingSlope * regionEtaMagnitude; ...
    -meanSwitchingSlope * regionEtaMagnitude; ...
    0.25 * ctrl.xi2 * regionEtaMagnitude];
for regionIndex = 0:3
    text(regionEta(regionIndex + 1), regionSigma(regionIndex + 1), ...
        sprintf('Region %d', regionIndex), 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 11, 'Color', [0.20 0.20 0.20], ...
        'BackgroundColor', [1.00 1.00 1.00], 'Margin', 2);
end
hold off; grid on; box on;
xlim([-etaExtent, etaExtent]);
xlabel('\eta = \int \sigma dt'); ylabel('\sigma = Sx');
title(sprintf('UAS-SMC phase plane (%s)', selectedLabel));
legend('Location', 'eastoutside');

figures.reachingLaw = figure('Name', 'UAS-SMC reaching law', 'Visible', visibility, ...
    'Color', 'w', 'Position', [180 180 1000 720]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(selected.t, selected.diagnostic.N, 'Color', [0.05 0.35 0.75], ...
    'LineWidth', 1.5); grid on; box on;
yline(0, '--k');
xlabel('Time (s)'); ylabel('N');
title(sprintf('Reaching law N(t) (%s)', selectedLabel));
nexttile;
plot(selected.t, selected.diagnostic.N_dot, 'Color', [0.80 0.20 0.15], ...
    'LineWidth', 1.3); grid on; box on;
yline(0, '--k');
xlabel('Time (s)'); ylabel('dN/dt');
title(sprintf('Numerical derivative of N(t), T_s = %.3f s', ...
    selected.t(2) - selected.t(1)));
end

function merged = localMergeStruct(defaults, overrides)
merged = defaults;
if isempty(overrides), return; end
names = fieldnames(overrides);
for idx = 1:numel(names)
    merged.(names{idx}) = overrides.(names{idx});
end
end
