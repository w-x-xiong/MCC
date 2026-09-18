function [x, fail, info] = SR_MCC(anc, rangems, maxiter, varargin)
% SR_MCC  Robust TOA localization via the squared-range maximum
% correntropy criterion (SR-MCC), solved by half-quadratic alternating
% maximization with a bisection-based GTRS inner solver.
%
% Reference: W. Xiong, C. Schindelhauer, H. C. So, Z. Wang, "Maximum
% Correntropy Criterion for Robust TOA-Based Localization in NLOS
% Environments," Circuits, Systems, and Signal Processing, 40,
% 6325-6339 (2021).
%
% This is a corrected and generalized re-implementation of the original
% MATLAB prototype (SR_MCC.m / bisection_fun.m). Differences are
% summarized below:
%
%   * Generalized from a hard-coded 2-D implementation (x = zeros(2,1);
%     D = zeros(3,3); ...) to arbitrary dimension H (2-D or 3-D).
%   * THE ORIGINAL sigma-UPDATE RULE NEVER ACTIVATES: the code
%     initializes sgm = inf and then sets
%         sgm = max(1.06*min(sgm_E,R/1.34)*L^(-0.2), sgm_old)
%     every iteration. Since max(anything finite, inf) = inf, sgm stays
%     at +inf for every iteration, all correntropy weights collapse to
%     one, and the algorithm silently reduces to *ordinary, non-robust*
%     squared-range least squares. This was confirmed empirically: the
%     original code's output is bit-identical (to solver tolerance) to
%     an unweighted SR-LS fit on every test case examined. Separately,
%     the kernel size sigma(0) = 0.1 stated in the paper's Algorithm-1
%     box is a fixed absolute constant while the squared-range residual
%     is in squared-length units; on typically-scaled problems this is
%     several orders of magnitude too small and collapses every weight
%     to numerically zero, producing NaN outputs (also confirmed
%     empirically). Both failure modes are fixed here by (a) a
%     data-driven, unit-consistent initial kernel size derived from an
%     unweighted warm-start fit, and (b) a schedule that is allowed to
%     shrink (subject to a floor and a maximum per-iteration shrink
%     rate), see schedule = 'annealed_mad' below.
%   * The GTRS solve is now dimension-general, uses a generalized
%     eigenvalue call instead of an explicit matrix square root, and
%     verifies/expands its bisection bracket (see wsrls_gtrs.m).
%
% Inputs
%   anc      : (H x L) sensor positions
%   rangems  : (1 x L) or (L x 1) TOA range measurements
%   maxiter  : maximum number of half-quadratic (AM) iterations
%
% Optional name-value pairs
%   'schedule'  : 'annealed_mad' (default, RECOMMENDED) |
%                 'fixed' | 'silverman_orig' | 'silverman_vanilla'
%                 ('silverman_orig' faithfully reproduces the original
%                  schedule's behaviour, sigma stuck at +inf, for
%                  side-by-side comparison; see note above.)
%   'sigmaFixed'  : value used when schedule = 'fixed'            (default 1)
%   'kappa0'      : initial-sigma multiple for 'annealed_mad'      (default 3)
%   'kappaMin'    : floor multiple for 'annealed_mad'              (default 0.1)
%   'eta'         : max per-iteration shrink factor, in (0,1)      (default 0.8)
%   'tol'         : convergence tolerance on ||x^(k+1)-x^(k)||     (default 1e-6)
%   'normalizeByRange' : if true, kernel argument is z_i/(2||x-x_i||)
%                        instead of z_i (see chapter remark on the
%                        heteroscedasticity of the SR residual)     (default false)
%
% Outputs
%   x     : (H x 1) location estimate
%   fail  : true if the algorithm failed to produce a finite estimate
%   info  : struct with fields iters, sigmaFinal, history (if requested)

p = inputParser;
addParameter(p, 'schedule', 'annealed_mad');
addParameter(p, 'sigmaFixed', 1.0);
addParameter(p, 'kappa0', 3.0);
addParameter(p, 'kappaMin', 0.1);
addParameter(p, 'eta', 0.8);
addParameter(p, 'tol', 1e-6);
addParameter(p, 'epsAbs', 1e-8);
addParameter(p, 'normalizeByRange', false);
addParameter(p, 'keepHistory', false);
parse(p, varargin{:});
opt = p.Results;

rangems = rangems(:)';
H = size(anc, 1);
L = size(anc, 2);

A1 = zeros(L, H + 1);
b1 = zeros(L, 1);
for i = 1:L
    A1(i, :) = [-2 * anc(:, i)', 1];
    b1(i) = rangems(i)^2 - norm(anc(:, i))^2;
end
D = zeros(H + 1, H + 1);
D(1:H, 1:H) = eye(H);
f = zeros(H + 1, 1);
f(end) = -0.5;

    function z = resid(xv)
        z = (rangems.^2 - sum((xv - anc).^2, 1))';
    end

    function a = effArg(xv, z)
        if opt.normalizeByRange
            dist = sqrt(sum((xv - anc).^2, 1))';
            a = z ./ max(2 * dist, 1e-6);
        else
            a = z;
        end
    end

% ---- warm start: unweighted SR-LS, then data-driven initial sigma ----
[x0, ~, ~, ok0] = wsrls_gtrs(A1, b1, D, f, ones(L, 1));
if ~ok0
    x0 = zeros(H, 1);
end
arg0 = effArg(x0, resid(x0));
scale0 = max(mad_scale(arg0), opt.epsAbs);

switch opt.schedule
    case 'fixed'
        sgm = opt.sigmaFixed;
    case {'silverman_orig', 'silverman_vanilla'}
        sgm = inf;
    case 'annealed_mad'
        sgm = max(opt.kappa0 * scale0, opt.epsAbs);
    otherwise
        error('SR_MCC:badSchedule', 'unknown schedule %s', opt.schedule);
end
sigma_min = max(opt.kappaMin * scale0, opt.epsAbs);

x = x0;
fail = false;
history = struct('k', {}, 'x', {}, 'sigma', {}, 'step', {});
k_final = 0;
for k = 1:maxiter
    k_final = k;
    x_old = x;

    z = resid(x);
    a = effArg(x, z);
    if ~isfinite(sgm) || sgm == 0
        w = ones(L, 1);
    else
        w = exp(-(a.^2) / (2 * sgm^2));
    end

    [x_new, ~, ~, ok] = wsrls_gtrs(A1, b1, D, f, w);
    if ~ok
        fail = true;
        break
    end

    z_new = resid(x_new);
    a_new = effArg(x_new, z_new);

    switch opt.schedule
        case 'fixed'
            % sgm unchanged
        case 'silverman_orig'
            % Faithful, literal reproduction of the original update rule
            % sgm = max(1.06*min(sgm_E,R/1.34)*L^(-0.2), sgm_old) with
            % sgm initialized to +inf. Because max(finite, Inf) = Inf in
            % IEEE arithmetic, sgm remains +inf for every iteration and
            % the correntropy weights never depart from 1 -- see the
            % file header for the empirically-confirmed consequence.
            s = silverman_scale(z_new, L);
            sgm = max(s, sgm);
        case 'silverman_vanilla'
            sgm = max(silverman_scale(a_new, L), opt.epsAbs);
        case 'annealed_mad'
            s_hat = mad_scale(a_new);
            target = min(sgm, max(opt.eta * sgm, s_hat));
            sgm = max(sigma_min, target);
    end

    step = norm(x_new - x_old);
    if opt.keepHistory
        history(end + 1) = struct('k', k, 'x', x_new, 'sigma', sgm, 'step', step); %#ok<AGROW>
    end

    x = x_new;
    if step < opt.tol
        break
    end
end

info.iters = k_final;
info.sigmaFinal = sgm;
info.history = history;

end

% -------------------------------------------------------------------
function s = mad_scale(z)
% 1.4826 * median absolute deviation: consistent estimator of sigma
% under Gaussian errors; 50% breakdown point.
s = 1.4826 * median(abs(z - median(z)));
end

function s = silverman_scale(z, L)
% Silverman (1986) rule-of-thumb bandwidth used in the original
% prototype: 1.06 * min(std, IQR/1.34) * L^{-1/5}. Implemented without
% the Statistics Toolbox's iqr() to avoid that dependency.
sd = std(z);
sz = sort(z(:));
q25 = prctile_manual(sz, 25);
q75 = prctile_manual(sz, 75);
iqrv = q75 - q25;
s = 1.06 * min(sd, iqrv / 1.34) * (L ^ (-0.2));
end

function q = prctile_manual(sorted_z, p)
% Linear-interpolation percentile (MATLAB/NumPy default-compatible),
% implemented manually so this file has no toolbox dependency.
n = numel(sorted_z);
if n == 1
    q = sorted_z(1);
    return
end
rank = (p / 100) * (n - 1) + 1;
lo = floor(rank);
hi = ceil(rank);
frac = rank - lo;
q = sorted_z(lo) + frac * (sorted_z(hi) - sorted_z(lo));
end
