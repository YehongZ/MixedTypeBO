function [model] = multigpMixedEP(model, iter)

% EP approximation for the non-Gaussian likelihood in multigpMixed
% Input:
%   model : the model structure
%	iter : maximal iterations of EP
% Output:
%   model : the model after doing EP approximation

    if numel(find(strcmp(model.oType, 'binary'))) == 0
        disp('No binary outputs for doing EP approximation.')
        return
    end

    persistent last_ttau last_tnu              % keep tilde parameters between calls
    tol = 1e-4; max_sweep = iter; min_sweep = 2;     % tolerance to stop EP iterations
   
    [K, mu0, y, ind, nd] = multigpMixedExtractBinary(model);
    n = length(mu0);

    
    [~, ~, ~, nlZ0] = epComputeParams(K, y, zeros(n, 1), zeros(n, 1), mu0);
    
    if any(size(last_ttau) ~= [n 1])      % find starting point for tilde parameters
        last_ttau = zeros(n, 1);
        last_tnu = zeros(n,1); 
    end

    ttau = last_ttau; tnu  = last_tnu;   % try the tilde values from previous call
    [Sigma,mu,~,nlZ] = epComputeParams(K, y, ttau, tnu, mu0);

    if nlZ > nlZ0                                           % if zero is better ..
        ttau = zeros(n,1); tnu  = zeros(n,1);       % .. then init with zero instead
        Sigma = K;                   % initialize Sigma and mu, the parameters of ..
        mu = mu0;
        nlZ = nlZ0;                % .. the Gaussian posterior approximation
    end

    nlZ_old = Inf; sweep = 0;               % converged, max. sweeps or min. sweeps?
    while (abs(nlZ-nlZ_old) > tol && sweep < max_sweep) || sweep<min_sweep

        nlZ_old = nlZ; sweep = sweep+1;

        for i = randperm(n)       % iterate EP updates (in random order) over examples
            
            tau_ni = 1/Sigma(i,i)-ttau(i);      %  first find the cavity distribution ..
            nu_ni = mu(i)/Sigma(i,i)-tnu(i);                % .. params tau_ni and nu_ni
            mu_ni = nu_ni/tau_ni;
                                     
            zi = y(i)*mu_ni/sqrt(1+1/tau_ni);
            ratio = exp((-0.5 * log(2 * pi) - 0.5 * zi.^2) - logphi(zi));

            mu_i = mu_ni + ratio*y(i)/sqrt(tau_ni*(1+tau_ni)); %%%%%%%%%%%%%%%
            sigma2i = 1/tau_ni - ratio/(tau_ni*(1+tau_ni))*(zi+ratio);

            dtt = 1/sigma2i - tau_ni - ttau(i);
            ttau_old = ttau(i);
            
            ttau(i) = ttau(i) + dtt;
            ttau(i) = max(ttau(i),0);
            tnu(i) = mu_i/sigma2i - nu_ni;
            
            dtt = ttau(i)-ttau_old;
            
            si = Sigma(:,i);
            ci = dtt/(1+dtt*si(i));
            Sigma = Sigma - ci*(si*si');   
            
            [Kinv, ~, ~] = pdinv(K);
            mu = Sigma*(tnu + Kinv*mu0);
        end
        
        [Sigma,mu,~,nlZ] = epComputeParams(K, y, ttau, tnu, mu0);
%         disp(['EP Step ', num2str(sweep), ': negative loglikelihood ', num2str(nlZ)])
    end

    if sweep == max_sweep && abs(nlZ-nlZ_old) > tol
        warning('maximum number of sweeps exceeded')
    end

    last_ttau = ttau; last_tnu = tnu;            % remember for next call
    
    ttau(ttau == 0) = 1;
    tsigma2 = 1./ttau;
    tmu = tnu./ttau; 
    N = 0;
    for i = 1:length(ind)
        s = length(model.oY{ind(i)+model.nlf});
        model.noise{ind(i)} = tsigma2(N+1:N+s);
        model.y(nd(i)+1:nd(i)+s) = tmu(N+1:N+s); 
        N = N + s;
        
    end   
    
    model = modelExpandParam(model, model.params);

% function to compute the parameters of the Gaussian approximation, Sigma and
% mu, and the negative log marginal likelihood, nlZ, from the current site
% parameters, ttau and tnu. Also returns L (useful for predictions).
function [Sigma,mu,L,nlZ] = epComputeParams(K, y, ttau, tnu, mu0)
    n = length(y);                                      % number of training cases
    sW = sqrt(ttau);
    L = chol(eye(n)+sW*sW'.*K);                            % L'*L=B=eye(n)+sW*K*sW
    V = L'\(repmat(sW,1,n).*K);
    Sigma = K - V'*V;
    
    [Kinv, ~, ~] = pdinv(K);
    mu = Sigma*(tnu + Kinv*mu0);   

    v = diag(Sigma);
    tau_n = 1./diag(Sigma)-ttau;             % compute the log marginal likelihood
    nu_n  = mu./diag(Sigma)-tnu;                    % vectors of cavity parameters
    lZ = logphi(y.*nu_n./tau_n./sqrt(1+1./tau_n));
    
    nlZ = sum(log(diag(L))) - sum(log(1+ttau./tau_n))/2 - sum(lZ) ...
        - tnu'*Sigma*tnu/2 + (v'*tnu.^2)/2 ...
            - nu_n'*((ttau./tau_n.*nu_n-2*tnu).*v)/2;
    
    V1 = L'\diag(sW);
    if any(ttau == 0)
        nlZ_mu0 = (mu0'/2)*(V1'*V1)*mu0;
    else
        nlZ_mu0 = -(tnu'./ttau' - mu0'/2)*(V1'*V1)*mu0;
    end
    nlZ = nlZ + nlZ_mu0;
    
