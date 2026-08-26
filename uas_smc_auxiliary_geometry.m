function geometry = uas_smc_auxiliary_geometry(ctrl)
%UAS_SMC_AUXILIARY_GEOMETRY Paper-consistent auxiliary-surface geometry.
%   The four auxiliary surfaces are centrally symmetric and their adjacent
%   intersections lie on s1 = 0 or s2 = 0, as constructed in Figure 2.6
%   and Equation (2.27) of the supplied UAS-SMC thesis.

required = {'xi1', 'xi2', 'beta'};
for idx = 1:numel(required)
    field = required{idx};
    if ~isfield(ctrl, field) || ~(isscalar(ctrl.(field)) && isfinite(ctrl.(field)))
        error('uas_smc:InvalidAuxiliaryGeometry', ...
            'Controller field %s must be a finite scalar.', field);
    end
end
if ctrl.xi1 <= ctrl.xi2 || ctrl.xi2 <= 0 || ctrl.beta <= 0
    error('uas_smc:InvalidAuxiliaryGeometry', ...
        'Require xi1 > xi2 > 0 and beta > 0.');
end
if isfield(ctrl, 'm_aux')
    mAux = ctrl.m_aux;
else
    mAux = 0.01;
end
if ~(isscalar(mAux) && isfinite(mAux) && mAux > 0)
    error('uas_smc:InvalidAuxiliaryGeometry', ...
        'm_aux must be a positive finite scalar.');
end

denominator = ctrl.xi1 + ctrl.xi2 + 2 * ctrl.beta;
alpha = (ctrl.beta * (ctrl.xi1 + ctrl.xi2) ...
    + 2 * ctrl.xi1 * ctrl.xi2) / denominator;
sideScale = (ctrl.xi1 - ctrl.xi2) / denominator;

% Rows correspond to regions 0, 1, 2, and 3.  The auxiliary surfaces are
% h_i = omega1(i)*sigma + omega2(i)*eta + m_aux.
omega1 = [1, -sideScale, sideScale, -1];
omega2 = [alpha, sideScale * ctrl.beta, ...
    -sideScale * ctrl.beta, -alpha];
coefficients = [omega1(:), omega2(:)];

% Adjacent faces meet in this counter-clockwise order.  Coordinates are
% returned as [sigma, eta].
facePairs = [2 4; 4 3; 3 1; 1 2];
vertices = zeros(4, 2);
for idx = 1:size(facePairs, 1)
    pair = facePairs(idx, :);
    vertices(idx, :) = (coefficients(pair, :) \ (-mAux * ones(2, 1))).';
end

s1 = vertices(:, 1) + ctrl.xi1 * vertices(:, 2);
s2 = vertices(:, 1) + ctrl.xi2 * vertices(:, 2);
switchingResidual = [s1(1); s2(2); s1(3); s2(4)];

geometry = struct( ...
    'alpha', alpha, ...
    'beta', ctrl.beta, ...
    'a', -alpha, ...
    'b', ctrl.beta, ...
    'sideScale', sideScale, ...
    'omega1', omega1, ...
    'omega2', omega2, ...
    'coefficients', coefficients, ...
    'm_aux', mAux, ...
    'vertices', vertices, ...
    'facePairs', facePairs, ...
    'switchingResidual', switchingResidual);
end
