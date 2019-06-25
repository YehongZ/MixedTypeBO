% Compute the recommended target maximizer (i.e., the maximizer of the
% predictive mean function)
function [f_xstar, xstar] = getMaxMean1(model, xmin, xmax, task, t, opts, x_prev)

    M = model.M;
    task_type = task{M+1};

    mu_func = @(x) getPosteriorMean(model, x, t);
    targetnew = @(x)(-mu_func(x));
    
    if opts.discrete
        Xtest = sample_sobol(size(xmin, 2), 1000, xmin, xmax);
        mu_f = getPosteriorMean(model, Xtest, t);
        [~, ind] = max(mu_f);
        xstart = Xtest(ind, :);
    else    
        bounds = [xmin', xmax'];
        problem.f = @(x)targetnew(x');
        opts.direct.maxevals = 1000;
        opts.cirect.maxits = 1000;

        [~, xstart] = Direct(problem, bounds, opts.direct);
        xstart = xstart';
    end
    
    [ x, fmin] = fmincon(targetnew, xstart, [], [], [], [], xmin, xmax, [], ...
                    optimset('MaxFunEvals', 500, 'AlwaysHonorConstraints', 'none', 'TolX', 1e-10, 'Tolfun', 1e-10, 'FinDiffType', 'central', 'Display', 'off'));
    
    % Compare with the recommended target maximizer from the previous BO
    % iteration, and select the one which has higher predictive mean.
    if size(size(x_prev), 2) > 2
        x_prev = reshape(x_prev, 1, size(x_prev, 3));
        x_prev = transX(x_prev, lower(task_type));

        mu_f = getPosteriorMean(model, x_prev, t);
        if mu_f > -fmin
            xstart = x_prev;
            [ x, ~] = fmincon(targetnew, xstart, [], [], [], [], xmin, xmax, [], ...
                    optimset('MaxFunEvals', 500, 'AlwaysHonorConstraints', 'none', 'TolX', 1e-10, 'Tolfun', 1e-10, 'FinDiffType', 'central', 'Display', 'off'));
        end
    end
    
    xstar = transX(x, lower(task_type), true);
    f_xstar = getFuncValue(task, xstar, M, t);
    
end

function [ ret ] = getPosteriorMean(model, xtest, t)

    nlf = model.nlf;
    [mu, ~] = multigpPosteriorMeanVar(model, xtest);
    ret = mu{t+nlf};

end


