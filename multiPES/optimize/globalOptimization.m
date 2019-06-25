function [ optimum, fval] = globalOptimization(target, gradient, xmin, xmax, opts)
    
    if opts.discrete
        num = opts.dis_num;        
        if size(xmin, 1) == 2
            a = linspace(xmin(1), xmax(1), num);
            b = linspace(xmin(2), xmax(2), num);
            [a1, b1] = meshgrid(a, b);
            Xtest = [reshape(a1, num*num, 1), reshape(b1, num*num, 1)];
        else
            Xtest = sample_sobol(size(xmin, 1), 1000, xmin', xmax');
        end
    end
    
    % Maximization -> Minimization
    targetnew = @(x) (-target(x));
    gradient = @(x) (-gradient(x));
    
    if opts.discrete
        value = target(Xtest');
        [~, index] = max(value);
        xstart = Xtest(index, :)';
    else
        % Direct - global optimization
        bounds = [xmin, xmax];
        problem.f = targetnew;
        opts.maxevals = 500;
        opts.maxits = 500;

        [~, xstart] = Direct(problem, bounds, opts.direct);
    end

	targetOptimization = @(x) deal(targetnew(x), gradient(x));

	% We optimize starting at the best location

	[ optimum, fval] = fmincon(targetOptimization, xstart, [], [], [], [], xmin, xmax, [], ...
                optimset('Algorithm', 'trust-region-reflective', 'MaxFunEvals', 500, 'Display', 'off', 'GradObj', 'on'));
           
end
