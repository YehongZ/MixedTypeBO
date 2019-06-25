function  ll = multigpMixedLogLikelihood(model)

% Compute the log likelihood of a multigpMixed model.
%
% COPYRIGHT : Mauricio A Alvarez, 2008
% This function is modified from MULTIGPLOGLIKELIHOOD by Yehong for fitting the
% mixed-type CMOGP algorithm.

[K, ~, y, ~, ~, ind] = multigpMixedExtractBinary_all(model);
ttau_all = 1./diag(model.K-model.fK);
tnu_all = ttau_all.*model.y;

n = length(ttau_all);                                      % number of training cases
sW = sqrt(ttau_all);
[L, ~] = chol(eye(n)+sW*sW'.*K);                            % L'*L=B=eye(n)+sW*K*sW
V = L'\(repmat(sW,1,n).*K);
Sigma = K - V'*V;
mu = Sigma*tnu_all;

v = diag(Sigma);
v = v(ind);
tau_n = 1./diag(Sigma)-ttau_all;             % compute the log marginal likelihood
nu_n  = mu./diag(Sigma)-tnu_all;
tau_n = tau_n(ind);
nu_n = nu_n(ind);
% vectors of cavity parameters
%p = tnu-m.*ttau; q = nu_n-m.*tau_n;                        % auxiliary vectors   
lZ = logphi(y.*nu_n./tau_n./sqrt(1+1./tau_n));

ttau = ttau_all(ind);
tnu = tnu_all(ind);

ttau_c = ttau_all(setdiff(1:length(ttau_all), ind));
tnu_c = tnu_all(setdiff(1:length(tnu_all), ind));

nlZ = sum(log(diag(L))) + sum(log(1./ttau_c))/2 - sum(log(1+ttau./tau_n))/2 - sum(lZ) ...
    - tnu_all'*Sigma*tnu_all/2 + (v'*tnu.^2)/2 + sum(tnu_c.^2./ttau_c/2) ...
        - nu_n'*((ttau./tau_n.*nu_n-2*tnu).*v)/2;

% nlZ = sum(log(diag(L))) - sum(log(1+ttau./tau_n))/2 - sum(lZ) ...
% - tnu'*Sigma*tnu/2 + (v'*tnu.^2)/2 ...
%     - nu_n'*((ttau./tau_n.*nu_n-2*tnu).*v)/2;

mu0 = model.y - model.m;
% mu0

V1 = L'\diag(sW);
if any(ttau_all == 0)
    nlZ_mu0 = (mu0'/2)*(V1'*V1)*mu0;
else
    nlZ_mu0 = -(tnu_all'./ttau_all' - mu0'/2)*(V1'*V1)*mu0;
end
nlZ = nlZ + nlZ_mu0;
    
ll = -nlZ;


    
