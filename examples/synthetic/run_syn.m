clear all, close all

Name = 'syn';
experimentNo = '_cr50';
cost = [50, 1];

N = 2500;   % budget
M = 2;  % no. of output types
d = 2;  % input dimension

T = 10; % repeat the experiments for 10 random initializations
opts.discrete = 1;
opts.dis_num = 10;
opts.direct.showits = 0;

xmin = [0, 0];
xmax = [1, 1];

t = 1;  % index of the target function
nSample = 50;
nFeatures = 200;
% bias m for generating the binary auxiliary function
me = [1.2, 1.1, 1.2, 1.2, 0.9, 0.8, 0.6, 1, 0.9, 0.5];  

% hypers of the synthetic functions
params = [log(100), log(100), log(1), log(2000), log(100), 1,...
log(100), log(2000), 1, log(0.01), log(0.001)];

for seq = 1:10
       
    task = {@target_synthetic, @auxiliary_synthetic, Name, seq, me(seq)};
    
    options = multigpOptions('ftc');
    options.oType = {'conti', 'binary'};
    options.M = M;
    options.q = d;
    options.params = params;
    options.nlf = 1;  % no. of latent function
    options.bias = [0, 0, -me(seq)];
    options.train = 0;
    
    result.X = zeros(T, N, d);
    result.y = zeros(T, N);
    result.type = zeros(T, N);
    result.cost = zeros(T, N);
    result.umax_f = zeros(T, N);
    result.umax_x = zeros(T, N, d);

    initX = cell(M, 1);
    
    for loop = 1:T
        
        rng('shuffle');
        
        for i=1:M
            initX{i} = rand(1, 2);
        end
       
        [model, Xobs, Yobs] = initializeBO(options, 'givenX', task, xmin, xmax, initX);
        S = sum(cost);
        
        [result.umax_f(loop, S-cost(2):S), result.umax_x(loop, 1, :)] = getMaxMean(model, xmin, xmax, task, t, opts);
        
        it = 1;
        while S <= N-cost(t)
            tic
            disp(['%%%%%%%%%%%%%%%%%%%%% Step: ', num2str(loop), ', ', num2str(it),  ', cost: ', num2str(S),' %%%%%%%%%%%%%%%%%%%%%']);
            [xstar, values, params] = sampleMixedXstar(model, nSample, nFeatures, xmin, xmax);
            
            disp('xstar sampling finished.');
            target_xstar = reshape(xstar(1, :, :), nSample, d);
            
            EPc = getFactor_EP(target_xstar, values, params, model.bias(model.nlf+1:end)); 

            [fi_xstar, modelnew] = EPapproximationMixed(target_xstar, model, 100, EPc); % Eq. (5.16)

            disp('First EP approximation finished.');
            acq_func = cell(M, 1);
            entropy_given_xstar = cell(M, 1);
            c = cell(M, 1);

            for i=1:M
                c{i} = getFactor(values, fi_xstar.m(:, i), i, t, 'avg', model, xstar);
                entropy_given_xstar{i} = @(x) EPapproximationMixed2(x, i, fi_xstar.m(:, i), fi_xstar.v(:, i), modelnew, c{i});    
                acq_func{i} = @(x) (entropy_cgp(model, x', i) - entropy_given_xstar{i}(x'))/cost(i);
            end    
           
            [optimum, type] = globalOptimizationNoGradient(acq_func, xmin', xmax', opts);
            b = toc;
            disp(['Next point selected: Elapsed time is ', num2str(b), ' seconds']);
          
            y = getObsValue(task, optimum', M, type); 
            [Xnew, ynew, noise] = updateMixedXY(model, optimum', y, type);

            tic
            model = updateBO_EP(model, Xnew, ynew, noise);
            b = toc;
            disp(['GP updated: Elapsed time is ', num2str(b), ' seconds']);
            
            model.xmin = xmin;
            model.xmax = xmax;
    
            result.X(loop, it, :) = optimum';
            result.y(loop, it) = y;
            result.type(loop, it) = type;
            result.cost(loop, it) = S + cost(type);
            result.umax_f(loop, S+1:S+cost(type)) = result.umax_f(loop, S);
            % umax_f is the true function value for input umax_x
            [result.umax_f(loop, S+cost(type)), result.umax_x(loop, it+1, :)] = getMaxMean1(model, xmin, xmax, task, t, opts, result.umax_x(loop, it, :));
            
            S = S + cost(type);
            it = it + 1;
        end 

        result.umax_f(loop, S+1:N) = result.umax_f(loop, S);
        save([Name, num2str(seq), experimentNo, '_result'], 'result', 'model');
    end
end
