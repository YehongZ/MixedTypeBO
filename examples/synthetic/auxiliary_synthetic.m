% Auxiliary synthetic function
% xx: 2-dimensional inputs in [0,1]^2
% seq: sequence number of the target function (1~10)
% m: bias of the function
% obs: 1 - return the binary output
%      0 - return the true function value
function [ ret ] = auxiliary_synthetic(xx, seq, m, obs)
    
    load(['synthetic', num2str(seq)]);
    
	[mu, ~] = multigpPosteriorMeanVar(model, xx);
    
    if nargin < 3
        m = 0;
    end
    
	ret = mu{3} - m;% + noise * randn(n, 1));
    
    if nargin < 4
        obs = false;
    end
    
    if obs
        ret = sign(exp(logphi(ret)) - 0.5);
    end
