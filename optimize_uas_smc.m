%OPTIMIZE_UAS_SMC Toolbox-free robust optimization of UAS-SMC parameters.
% Optional caller variable: optimizer_options (structure).  The script
% leaves optimization_report and best_ctrl in the caller workspace.

scriptRoot = fileparts(mfilename('fullpath'));
defaultOptions = struct( ...
    'particles', 40, ...
    'iterations', 60, ...
    'trainingRandomCount', 12, ...
    'validationCount', 100, ...
    'localMaxRounds', 100, ...
    'saveResults', true, ...
    'verbose', true, ...
    'outputDir', fullfile(scriptRoot, 'results'), ...
    'seed', 20260717, ...
    'validationSeed', 20260718, ...
    'Ts', 0.005, ...
    'Tend', 8, ...
    'inertiaWeight', 0.72, ...
    'cognitiveWeight', 1.49, ...
    'socialWeight', 1.49, ...
    'velocityFraction', 0.20);

if ~exist('optimizer_options', 'var') || isempty(optimizer_options)
    optimizer_options = defaultOptions;
else
    optimizer_options = localMergeStruct(defaultOptions, optimizer_options);
end

parameterNames = {'p1', 'poleGap12', 'poleGap23', 'xi2', 'xiGap', 'beta', 'k'};
lowerBounds = [0.5, 0.25, 0.25, 0.2, 0.2, 0.1, 0.2];
upperBounds = [5.0, 5.00, 5.00, 8.0, 12., 8.0, 15.];
baselineQ = [2, 1, 1, 1, 1, 1, 4];

nominalParams = localNominalParameters();
nominalModel = localBuildModel(nominalParams);
trainingScenarios = localCreateTrainingScenarios(nominalParams, optimizer_options);
trainingScenarios = localPrepareScenarios(trainingScenarios, optimizer_options.Ts);

if optimizer_options.verbose
    fprintf('UAS-SMC robust optimization\n');
    fprintf('  Training scenarios : %d\n', numel(trainingScenarios));
    fprintf('  Particles/iterations: %d / %d\n', ...
        optimizer_options.particles, optimizer_options.iterations);
end

baselineCtrl = localDecodeController(baselineQ);
[baselineScore, baselineEvaluation] = localObjective( ...
    baselineQ, nominalModel, trainingScenarios, optimizer_options);

rng(optimizer_options.seed, 'twister');
dimension = numel(lowerBounds);
range = upperBounds - lowerBounds;
particleCount = optimizer_options.particles;
positions = lowerBounds + rand(particleCount, dimension) .* range;
positions(1, :) = baselineQ;
velocityLimit = optimizer_options.velocityFraction * range;
velocities = (2 * rand(particleCount, dimension) - 1) .* velocityLimit;
velocities(1, :) = 0;

scores = inf(particleCount, 1);
evaluations = cell(particleCount, 1);
for particle = 1:particleCount
    [scores(particle), evaluations{particle}] = localObjective( ...
        positions(particle, :), nominalModel, trainingScenarios, optimizer_options);
end
personalPositions = positions;
personalScores = scores;
[globalScore, bestIndex] = min(scores);
globalPosition = positions(bestIndex, :);
globalEvaluation = evaluations{bestIndex};

historyBest = zeros(optimizer_options.iterations + 1, 1);
historyMean = zeros(optimizer_options.iterations + 1, 1);
historyBest(1) = globalScore;
historyMean(1) = mean(scores);

for iteration = 1:optimizer_options.iterations
    randomPersonal = rand(particleCount, dimension);
    randomGlobal = rand(particleCount, dimension);
    velocities = optimizer_options.inertiaWeight * velocities ...
        + optimizer_options.cognitiveWeight * randomPersonal .* (personalPositions - positions) ...
        + optimizer_options.socialWeight * randomGlobal .* (globalPosition - positions);
    velocities = min(max(velocities, -velocityLimit), velocityLimit);
    positions = min(max(positions + velocities, lowerBounds), upperBounds);

    for particle = 1:particleCount
        [scores(particle), evaluations{particle}] = localObjective( ...
            positions(particle, :), nominalModel, trainingScenarios, optimizer_options);
        if scores(particle) < personalScores(particle)
            personalScores(particle) = scores(particle);
            personalPositions(particle, :) = positions(particle, :);
        end
        if scores(particle) < globalScore
            globalScore = scores(particle);
            globalPosition = positions(particle, :);
            globalEvaluation = evaluations{particle};
        end
    end
    historyBest(iteration + 1) = globalScore;
    historyMean(iteration + 1) = mean(scores);
    if optimizer_options.verbose && (iteration == 1 || mod(iteration, 5) == 0)
        fprintf('  PSO iteration %3d/%3d: best=%.8g, mean=%.8g\n', ...
            iteration, optimizer_options.iterations, globalScore, historyMean(iteration + 1));
    end
end

step = 0.10 * range;
minimumStep = 0.001 * range;
localScores = globalScore;
localRounds = 0;
for roundIndex = 1:optimizer_options.localMaxRounds
    localRounds = roundIndex;
    candidateBestScore = globalScore;
    candidateBestPosition = globalPosition;
    candidateBestEvaluation = globalEvaluation;
    for variable = 1:dimension
        for direction = [-1, 1]
            candidate = globalPosition;
            candidate(variable) = min(max(candidate(variable) ...
                + direction * step(variable), lowerBounds(variable)), upperBounds(variable));
            if candidate(variable) == globalPosition(variable)
                continue;
            end
            [candidateScore, candidateEvaluation] = localObjective( ...
                candidate, nominalModel, trainingScenarios, optimizer_options);
            if candidateScore < candidateBestScore
                candidateBestScore = candidateScore;
                candidateBestPosition = candidate;
                candidateBestEvaluation = candidateEvaluation;
            end
        end
    end
    if candidateBestScore < globalScore
        globalScore = candidateBestScore;
        globalPosition = candidateBestPosition;
        globalEvaluation = candidateBestEvaluation;
    else
        step = 0.5 * step;
    end
    localScores(end + 1, 1) = globalScore; %#ok<SAGROW>
    if all(step <= minimumStep)
        break;
    end
end

% The baseline is explicitly part of the candidate set, so the saved result
% can never be worse than it under the declared objective.
if baselineScore < globalScore
    globalScore = baselineScore;
    globalPosition = baselineQ;
    globalEvaluation = baselineEvaluation;
end

best_ctrl = localDecodeController(globalPosition);
validationScenarios = localCreateValidationScenarios(nominalParams, optimizer_options);
validationScenarios = localPrepareScenarios(validationScenarios, optimizer_options.Ts);
[validationScore, validationEvaluation] = localObjective( ...
    globalPosition, nominalModel, validationScenarios, optimizer_options);

optimization_report = struct();
optimization_report.metadata = struct( ...
    'seed', optimizer_options.seed, ...
    'validationSeed', optimizer_options.validationSeed, ...
    'particles', optimizer_options.particles, ...
    'iterations', optimizer_options.iterations, ...
    'trainingScenarioCount', numel(trainingScenarios), ...
    'validationScenarioCount', numel(validationScenarios), ...
    'parameterCount', numel(parameterNames), ...
    'parameterNames', {parameterNames}, ...
    'localRounds', localRounds, ...
    'Ts', optimizer_options.Ts, ...
    'Tend', optimizer_options.Tend);
optimization_report.bounds = struct('lower', lowerBounds, 'upper', upperBounds);
optimization_report.baseline = struct( ...
    'q', baselineQ, 'ctrl', baselineCtrl, 'score', baselineScore, ...
    'evaluation', baselineEvaluation);
optimization_report.best = struct( ...
    'q', globalPosition, 'ctrl', best_ctrl, 'score', globalScore, ...
    'evaluation', globalEvaluation);
optimization_report.validation = struct( ...
    'score', validationScore, ...
    'evaluation', validationEvaluation, ...
    'passRate', mean(validationEvaluation.feasible));
optimization_report.history = struct( ...
    'bestScore', historyBest, ...
    'meanScore', historyMean, ...
    'localBestScore', localScores);

if optimizer_options.saveResults
    if exist(optimizer_options.outputDir, 'dir') ~= 7
        mkdir(optimizer_options.outputDir);
    end
    best_q = globalPosition;
    save(fullfile(optimizer_options.outputDir, 'uas_smc_best_params.mat'), ...
        'best_ctrl', 'best_q', 'optimization_report');
    save(fullfile(optimizer_options.outputDir, 'optimization_history.mat'), ...
        'optimization_report');
    optimizationFigure = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 900 560]);
    semilogy(0:optimizer_options.iterations, historyBest, 'LineWidth', 1.6);
    hold on;
    semilogy(0:optimizer_options.iterations, historyMean, '--', 'LineWidth', 1.1);
    if numel(localScores) > 1
        localAxis = optimizer_options.iterations + (0:numel(localScores)-1);
        semilogy(localAxis, localScores, 'LineWidth', 1.4);
    end
    grid on;
    xlabel('Optimization iteration');
    ylabel('Robust objective');
    legend('PSO best', 'PSO population mean', 'Local refinement', 'Location', 'best');
    title('UAS-SMC optimization convergence');
    exportgraphics(optimizationFigure, ...
        fullfile(optimizer_options.outputDir, 'optimization_convergence.png'), ...
        'Resolution', 180);
    close(optimizationFigure);
end

if optimizer_options.verbose
    fprintf('\nOptimization result\n');
    fprintf('  Baseline / best score : %.8g / %.8g\n', baselineScore, globalScore);
    fprintf('  Best q                : [%s]\n', num2str(globalPosition, ' %.7g'));
    fprintf('  Internal poles        : [%s]\n', num2str(-best_ctrl.poles, ' %.7g'));
    fprintf('  xi1, xi2, a, b, k     : %.7g, %.7g, %.7g, %.7g, %.7g\n', ...
        best_ctrl.xi1, best_ctrl.xi2, best_ctrl.a, best_ctrl.b, best_ctrl.k);
    fprintf('  Training pass rate    : %.1f%%\n', 100 * mean(globalEvaluation.feasible));
    fprintf('  Validation pass rate  : %.1f%%\n\n', ...
        100 * mean(validationEvaluation.feasible));
end

function params = localNominalParameters()
params = struct('m', 0.0923, 'M', 0.1945, 'l', 0.1950, 'J', 0.0029, ...
    'g', 9.8, 'R', 3.75, 'r', 0.018, 'Km', 0.183, 'Ke', 0.232, 'I', 7.083e-6);
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
model = struct('params', p, 'A', A, 'B', B, 'Co', Co);
end

function ctrl = localDecodeController(q)
p1 = q(1);
p2 = p1 + q(2);
p3 = p2 + q(3);
ctrl = struct( ...
    'poles', [p1 p2 p3], ...
    'xi2', q(4), ...
    'xi1', q(4) + q(5), ...
    'beta', q(6), ...
    'k', q(7), ...
    'm_aux', 0.01, ...
    'u_max', 12, ...
    'tol', 1e-12);
geometry = uas_smc_auxiliary_geometry(ctrl);
ctrl.alpha = geometry.alpha;
ctrl.a = geometry.a;
ctrl.b = geometry.b;
ctrl.auxiliary_side_scale = geometry.sideScale;
end

function design = localBuildDesign(model, ctrl)
h = [0 0 0 1] / model.Co;
transform = [h; h * model.A; h * model.A^2; h * model.A^3];
polynomial = poly(-ctrl.poles);
surfaceCoefficients = [polynomial(4), polynomial(3), polynomial(2), 1];
S = surfaceCoefficients * transform;
S = S / (S * model.B);
design = struct('S', S, 'SB', S * model.B);
end

function scenarios = localCreateTrainingScenarios(nominal, options)
fixedAngles = [5, -5, 8, -8];
scenarios = repmat(struct('params', nominal, 'x0', zeros(4, 1), 'label', ''), ...
    1, numel(fixedAngles) + options.trainingRandomCount);
for idx = 1:numel(fixedAngles)
    scenarios(idx).params = nominal;
    scenarios(idx).x0 = [0; 0; deg2rad(fixedAngles(idx)); 0];
    scenarios(idx).label = sprintf('nominal_%+gdeg', fixedAngles(idx));
end
rng(options.seed, 'twister');
for idx = 1:options.trainingRandomCount
    scenarioIndex = numel(fixedAngles) + idx;
    scenarios(scenarioIndex).params = localPerturbParameters(nominal);
    angle = -8 + 16 * rand();
    scenarios(scenarioIndex).x0 = [0; 0; deg2rad(angle); 0];
    scenarios(scenarioIndex).label = sprintf('training_random_%02d', idx);
end
end

function scenarios = localCreateValidationScenarios(nominal, options)
scenarios = repmat(struct('params', nominal, 'x0', zeros(4, 1), 'label', ''), ...
    1, options.validationCount);
rng(options.validationSeed, 'twister');
for idx = 1:options.validationCount
    scenarios(idx).params = localPerturbParameters(nominal);
    angle = -8 + 16 * rand();
    scenarios(idx).x0 = [0; 0; deg2rad(angle); 0];
    scenarios(idx).label = sprintf('validation_%03d', idx);
end
end

function perturbed = localPerturbParameters(nominal)
perturbed = nominal;
fields = {'M', 'm', 'J', 'l', 'R', 'Ke', 'Km', 'I'};
for idx = 1:numel(fields)
    perturbed.(fields{idx}) = nominal.(fields{idx}) * (0.9 + 0.2 * rand());
end
end

function scenarios = localPrepareScenarios(scenarios, Ts)
for idx = 1:numel(scenarios)
    plant = localBuildModel(scenarios(idx).params);
    augmentedA = [plant.A, zeros(4); eye(4), zeros(4)];
    augmentedB = [plant.B; zeros(4, 1)];
    transition = expm([augmentedA, augmentedB; zeros(1, 9)] * Ts);
    scenarios(idx).A = plant.A;
    scenarios(idx).B = plant.B;
    scenarios(idx).Ad = transition(1:4, 1:4);
    scenarios(idx).Bd = transition(1:4, 9);
    scenarios(idx).integralStateGain = transition(5:8, 1:4);
    scenarios(idx).integralInputGain = transition(5:8, 9);
end
end

function [score, evaluation] = localObjective(q, nominalModel, scenarios, options)
ctrl = localDecodeController(q);
design = localBuildDesign(nominalModel, ctrl);
scenarioCount = numel(scenarios);
sampleCount = round(options.Tend / options.Ts) + 1;
state = zeros(4, scenarioCount);
eta = zeros(1, scenarioCount);
for scenario = 1:scenarioCount
    state(:, scenario) = scenarios(scenario).x0;
end

iaePosition = zeros(1, scenarioCount);
iaeAngle = zeros(1, scenarioCount);
sumVoltageSquared = zeros(1, scenarioCount);
sumDeltaVoltageSquared = zeros(1, scenarioCount);
peakPosition = zeros(1, scenarioCount);
peakAngle = zeros(1, scenarioCount);
peakVoltage = zeros(1, scenarioCount);
peakRawVoltage = zeros(1, scenarioCount);
minimumN = inf(1, scenarioCount);
previousVoltage = zeros(1, scenarioCount);
lastOutsideTime = -options.Ts * ones(1, scenarioCount);
allFinite = true(1, scenarioCount);
geometry = uas_smc_auxiliary_geometry(ctrl);

for sample = 1:sampleCount
    currentTime = (sample - 1) * options.Ts;
    sigma = design.S * state;
    s1 = sigma + ctrl.xi1 * eta;
    s2 = sigma + ctrl.xi2 * eta;

    omega1 = geometry.omega1(1) * ones(1, scenarioCount);
    omega2 = geometry.omega2(1) * ones(1, scenarioCount);
    region1 = s1 < 0 & s2 >= 0;
    region2 = s1 >= 0 & s2 < 0;
    region3 = s1 >= 0 & s2 >= 0;
    omega1(region1) = geometry.omega1(2);
    omega2(region1) = geometry.omega2(2);
    omega1(region2) = geometry.omega1(3);
    omega2(region2) = geometry.omega2(3);
    omega1(region3) = geometry.omega1(4);
    omega2(region3) = geometry.omega2(4);

    epsilon = ones(1, scenarioCount);
    firstBlend = s1 .* s2 <= 0 & abs(s1) > ctrl.tol;
    epsilon(firstBlend) = abs(s2(firstBlend)) ./ ...
        (abs(s1(firstBlend)) + abs(s2(firstBlend)));
    secondBlend = ~firstBlend & s2 .* sigma <= 0 & abs(sigma) > ctrl.tol;
    epsilon(secondBlend) = abs(s2(secondBlend)) ./ ...
        (abs(s2(secondBlend)) + abs(sigma(secondBlend)));

    phi = epsilon .* (ctrl.a * sigma - ctrl.k * s2) ...
        + (1 - epsilon) .* 0.5 * (ctrl.a + ctrl.b) .* sigma;
    reachingLaw = omega2 .* sigma + omega1 .* phi;
    rawVoltage = (-design.S * nominalModel.A * state + phi) / design.SB;
    voltage = min(max(rawVoltage, -ctrl.u_max), ctrl.u_max);

    position = state(1, :);
    angle = state(3, :);
    iaePosition = iaePosition + abs(position) * options.Ts;
    iaeAngle = iaeAngle + abs(angle) * options.Ts;
    sumVoltageSquared = sumVoltageSquared + voltage.^2;
    if sample > 1
        sumDeltaVoltageSquared = sumDeltaVoltageSquared + (voltage - previousVoltage).^2;
    end
    previousVoltage = voltage;
    peakPosition = max(peakPosition, abs(position));
    peakAngle = max(peakAngle, abs(angle));
    peakVoltage = max(peakVoltage, abs(voltage));
    peakRawVoltage = max(peakRawVoltage, abs(rawVoltage));
    minimumN = min(minimumN, reachingLaw);
    within = abs(position) < 1e-3 & abs(angle) < deg2rad(0.1);
    lastOutsideTime(~within) = currentTime;
    allFinite = allFinite & all(isfinite(state), 1) & isfinite(voltage) & isfinite(reachingLaw);

    if sample < sampleCount
        for scenario = 1:scenarioCount
            oldState = state(:, scenario);
            state(:, scenario) = scenarios(scenario).Ad * oldState ...
                + scenarios(scenario).Bd * voltage(scenario);
            stateIntegral = scenarios(scenario).integralStateGain * oldState ...
                + scenarios(scenario).integralInputGain * voltage(scenario);
            eta(scenario) = eta(scenario) + design.S * stateIntegral;
        end
    end
end

finalWithin = abs(state(1, :)) < 1e-3 & abs(state(3, :)) < deg2rad(0.1);
settlingTime = max(0, lastOutsideTime + options.Ts);
settlingTime(~finalWithin) = options.Tend;
rmsVoltage = sqrt(sumVoltageSquared / sampleCount);
rmsDeltaVoltage = sqrt(sumDeltaVoltageSquared / max(1, sampleCount - 1));

caseObjective = 0.30 * (settlingTime / options.Tend) ...
    + 0.15 * (iaePosition / (0.2 * options.Tend)) ...
    + 0.20 * (iaeAngle / (deg2rad(10) * options.Tend)) ...
    + 0.10 * (peakPosition / 0.2) ...
    + 0.10 * (rmsVoltage / 12) ...
    + 0.10 * (rmsDeltaVoltage / 12) ...
    + 0.05 * (peakAngle / deg2rad(10));

feasible = allFinite ...
    & peakAngle <= deg2rad(10) + 1e-12 ...
    & peakPosition <= 0.20 + 1e-12 ...
    & peakRawVoltage <= 12 + 1e-9 ...
    & finalWithin ...
    & minimumN >= -1e-9;

angleViolation = max(0, peakAngle / deg2rad(10) - 1);
positionViolation = max(0, peakPosition / 0.2 - 1);
voltageViolation = max(0, peakRawVoltage / 12 - 1);
finalViolation = max(0, abs(state(1, :)) / 1e-3 - 1) ...
    + max(0, abs(state(3, :)) / deg2rad(0.1) - 1);
reachingViolation = max(0, -minimumN / 1e-9);
penalty = 1e6 * sum(~feasible) + 1e4 * sum(angleViolation ...
    + positionViolation + voltageViolation + finalViolation + reachingViolation);
score = 0.6 * max(caseObjective) + 0.4 * mean(caseObjective) + penalty;
if ~isfinite(score)
    score = realmax('double') / 1e100;
end

evaluation = struct( ...
    'caseObjective', caseObjective, ...
    'feasible', feasible, ...
    'settlingTime', settlingTime, ...
    'iaePosition', iaePosition, ...
    'iaeAngle', iaeAngle, ...
    'peakPosition', peakPosition, ...
    'peakAngleDeg', rad2deg(peakAngle), ...
    'peakVoltage', peakVoltage, ...
    'peakRawVoltage', peakRawVoltage, ...
    'rmsVoltage', rmsVoltage, ...
    'rmsDeltaVoltage', rmsDeltaVoltage, ...
    'minimumN', minimumN, ...
    'allFinite', allFinite, ...
    'inputGain', design.SB, ...
    'finalPosition', state(1, :), ...
    'finalAngleDeg', rad2deg(state(3, :)));
end

function merged = localMergeStruct(defaults, overrides)
merged = defaults;
if isempty(overrides), return; end
names = fieldnames(overrides);
for idx = 1:numel(names)
    merged.(names{idx}) = overrides.(names{idx});
end
end
