function [ optimum, type, xs ] = globalOptimizationNoGradient(target, xmin, xmax, opts)
    
    num = opts.dis_num; 
    if opts.discrete           
        if size(xmin, 1) == 2
            a = linspace(xmin(1), xmax(1), num);
            b = linspace(xmin(2), xmax(2), num);
            [a1, b1] = meshgrid(a, b);
            Xtest = [reshape(a1, num*num, 1), reshape(b1, num*num, 1)];
        else
            Xtest = sample_sobol(size(xmin, 1), opts.num, xmin', xmax');
        end
    end
    
    fval = inf;
    for i=1:size(target, 1)
        % Maximization -> Minimization
        targetnew = @(x)(-target{i}(x));
        
        if opts.discrete
            tmp = -Inf;
            xstart = xmin;
            for j = 1:size(Xtest, 1) 
                value = target{i}(Xtest(j, :)');
                if value > tmp
                    tmp = value;
                    xstart = Xtest(j, :)';
                end
            end
        else
            % Direct - global optimization
            bounds = [xmin, xmax];
            problem.f = targetnew;
            opts.direct.maxevals = 500;
            opts.cirect.maxits = 500;

            [tmp, xstart] = Direct(problem, bounds, opts.direct);
            tmp = -tmp;
        end

        % We optimize starting at the best location
        [ x, fmin] = fmincon(targetnew, xstart, [], [], [], [], xmin, xmax, [], ...
                    optimset('MaxFunEvals', 500, 'AlwaysHonorConstraints', 'none', 'TolX', 1e-10, 'Tolfun', 1e-10, 'FinDiffType', 'central', 'Display', 'off'));
        
        if tmp > -fmin
            
            disp('----------------------')
            inval = (xmax-xmin)./num;
            xmin_new = max(xstart-inval, xmin);
            xmax_new = min(xstart+inval, xmax);
            [ x, fmin] = fmincon(targetnew, xstart, [], [], [], [], xmin_new, xmax_new, [], ...
                optimset('MaxFunEvals', 500, 'AlwaysHonorConstraints', 'none', 'TolX', 1e-10, 'Tolfun', 1e-10, 'FinDiffType', 'central', 'Display', 'off'));
            
        end
        disp(['Type: ', num2str(i), ' f_start: ', num2str(tmp), ' Acq_max: ' num2str(-fmin)])
        if fmin < fval
            optimum = x;
            fval = fmin;
            type = i;
            xs = xstart;         
        end 
    end
           
end

% 	optimum = fmincon(target, start, [], [], [], [], xmin, xmax, [], ...
%                 optimset('MaxFunEvals', 100, 'TolX', eps, 'Display', 'off', 'GradObj', 'on'));
