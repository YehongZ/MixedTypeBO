% Compute the slack variables c_i
function ret = getFactor_EP(xt, fxstar, params, bias)

    [nSample, M] = size(fxstar);
    astar = zeros(nSample, M);
    for i=1:nSample

        theta = params{i}.theta;
        trans_all = params{i}.trans;
        W_all = params{i}.W;
        b_all = params{i}.b;
        x = xt(i, :)';
        phi = cos(W_all * x + repmat(b_all, 1, size(x, 2)));
        theta = repmat(theta, 1, M) .* trans_all;
        
        astar(i, :) = phi'*theta + bias;

    end
    ret = mean(max(fxstar - astar, 0), 1)';
    
    