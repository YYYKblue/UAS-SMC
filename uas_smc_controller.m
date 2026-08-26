function [u, diagnostic] = uas_smc_controller(x, eta, model, ctrl)
%UAS_SMC_CONTROLLER Continuous UAS-SMC law for a scalar virtual state.
%   The controllable four-state plant is reduced to sigma = S*x.  The
%   switching surfaces, active auxiliary surface, and reaching law follow
%   the symmetric no-chattering construction in Equation (2.57) of the
%   supplied UAS-SMC paper.

x = x(:);
if numel(x) ~= 4
    error('uas_smc:InvalidState', 'x must contain four states.');
end
if ~(isscalar(eta) && isfinite(eta))
    error('uas_smc:InvalidIntegralState', 'eta must be a finite scalar.');
end
if ctrl.xi1 <= ctrl.xi2 || ctrl.xi2 <= 0
    error('uas_smc:InvalidSwitchingSurfaces', 'Require xi1 > xi2 > 0.');
end
if ctrl.beta <= 0 || ctrl.k <= 0
    error('uas_smc:InvalidControllerParameters', ...
        'Require beta and k to be positive.');
end
if abs(model.SB) <= ctrl.tol
    error('uas_smc:SingularInputChannel', 'The scalar input gain S*B is too small.');
end

sigma = model.S * x;
s1 = sigma + ctrl.xi1 * eta;
s2 = sigma + ctrl.xi2 * eta;
geometry = uas_smc_auxiliary_geometry(ctrl);

if s1 < 0 && s2 < 0
    region = 0;
elseif s1 < 0 && s2 >= 0
    region = 1;
elseif s1 >= 0 && s2 < 0
    region = 2;
else
    region = 3;
end
omega1 = geometry.omega1(region + 1);
omega2 = geometry.omega2(region + 1);

if s1 * s2 <= 0 && abs(s1) > ctrl.tol
    epsilon = abs(s2) / (abs(s1) + abs(s2));
elseif s2 * sigma <= 0 && abs(sigma) > ctrl.tol
    epsilon = abs(s2) / (abs(s2) + abs(sigma));
else
    epsilon = 1;
end

a = geometry.a;
b = geometry.b;
phi = epsilon * (a * sigma - ctrl.k * s2) ...
    + (1 - epsilon) * 0.5 * (a + b) * sigma;
N = omega2 * sigma + omega1 * phi;
h = omega1 * sigma + omega2 * eta + geometry.m_aux;

uUnsaturated = (-model.S * model.A * x + phi) / model.SB;
u = min(max(uUnsaturated, -ctrl.u_max), ctrl.u_max);

if any(~isfinite([u, sigma, s1, s2, epsilon, phi, N, h]))
    error('uas_smc:NonFiniteControl', 'The UAS-SMC calculation produced a non-finite value.');
end

diagnostic = struct( ...
    'sigma', sigma, ...
    's1', s1, ...
    's2', s2, ...
    'epsilon', epsilon, ...
    'phi', phi, ...
    'N', N, ...
    'h', h, ...
    'region', region, ...
    'omega1', omega1, ...
    'omega2', omega2, ...
    'u_unsat', uUnsaturated, ...
    'saturated', abs(uUnsaturated) > ctrl.u_max + ctrl.tol);
end
