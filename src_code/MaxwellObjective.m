function L2 = MaxwellObjective(md, params, t_maxwell)
% Objective function for inverting for input parameters for a desired
% Maxwell rheology

% compute desired relaxation functions of Maxwell material----------------
% unpack parameters
n = params(1); rho = params(2); eta = params(3);
mu = params(4); K = params(5); vp_0 = params(6);
vs_0 = params(7);

G1_M = 2*mu*exp(-t_maxwell/(eta/mu));
G2_M = repelem(3*K, length(t_maxwell), 1);

% compute approximate relaxation functions from SeisSol-------------------
% unpack variales
f_c = 10.^md(1); f_ratio = 10.^md(2); 
QP = 10.^md(3); QS = 10.^md(4);
vp_t = 10.^md(5); vs_t = 10.^md(6);

% first run to get the optimal vp_0, vs_0 (velocities at infinite frequency)
[~, ~, params] = cmp_relaxationFunc(QP, QS, n, f_c, f_ratio, vp_0, vs_0, vp_t, vs_t, rho);

% second run to get the corresponding relaxation functions
[G1, G2, ~] = cmp_relaxationFunc(QP, QS, n, f_c, f_ratio, params.vp_0, params.vs_0, vp_t, vs_t, rho);

G1_approx = G1(t_maxwell')'; % note input time has to be N x 1
G2_approx = G2(t_maxwell')';

L2 = norm(log(G1_M) - log(G1_approx)) + norm(log(G2_M) - log(G2_approx)); % log form is necessary to find good fit

end