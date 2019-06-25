% Target synthetic function
% xx: 2-dimensional inputs in [0,1]^2
% seq: sequence number of the target function (1~10)
% m: bias of the function
% obs: 1 - return the noisy output
%      0 - return the true function value
function [ ret ] = target_synthetic(xx, seq, m, obs)
    
    load(['synthetic', num2str(seq)]);
    noise = sqrt(exp(model.params(10)));
    
	[mu, ~] = multigpPosteriorMeanVar(model, xx);
    
	ret = mu{2};% + noise * randn(n, 1));
    
    if nargin < 4
        obs = false;
    end
    
    if obs
        n = size(xx, 1);
        ret = ret + noise * randn(n, 1);
    end
