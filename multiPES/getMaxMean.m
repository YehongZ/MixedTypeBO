% Compute the recommended target maximizer (i.e., the maximizer of the
% predictive mean function)
function [f_xstar, xstar] = getMaxMean(model, xmin, xmax, task, t, opts)
    
    if strcmp(model.oType{t}, 'binary') == 1
        error('Invalid output type for the primary task');
    end
    
    M = model.M;
    
    mu_func{1} = @(x) getPosteriorMean(model, x', t);
    
    opts.dis_num = 20;
    opts.discrete = 1;
    
    [xstar, ~] = globalOptimizationNoGradient(mu_func, xmin', xmax', opts);
    xstar = xstar';
    
    task_type = task{M+1};
    xstar = transX(xstar, lower(task_type), true);
    f_xstar = getFuncValue(task, xstar, M, t);
end

function [ ret ] = getPosteriorMean(model, xtest, t)

    nlf = model.nlf;
    [mu, ~] = multigpMixedPosteriorMeanVar(model, xtest);
    ret = mu{t+nlf};

end


