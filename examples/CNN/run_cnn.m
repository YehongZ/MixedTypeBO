clear all, close all

Name = 'CNN';
experimentNo = ['_bos05'];
fname = ['cnn', experimentNo, '.txt'];  % file used to transfer the CNN training results

% prior trained hypers for mixed-type CMOGP model
load('cnn_params');

gpu = 0;
N = 1000;
Budget = 36; % training budget: 36hr
M = 2;
d = 6;  % input dimension
T = 5;  % repeat the experiments for 5 random initializations

opts.discrete = 1;
opts.dis_num = 10;
opts.num = 500;
opts.direct.showits = 0;

task = {{}, {}, Name};

% Input: epochs, three dropout rate, learning rate, batch_size
xmin = [100, 0, 0, 0, log(10.^-5), 100];
xmax = [1000, 1, 1, 1, log(1), 1000];

t = 1;
nSample = 50;
nFeatures = 200;

result.X = zeros(T, N, d);
result.y = zeros(T, N);
result.type = zeros(T, N);
result.umax_f = zeros(T, N);
result.umax_x = zeros(T, N, d);
result.time_real = zeros(T, N);

load('start_cnn_bos');
initX = cell(M, 1);
initY = cell(M, 1);
initT = cell(M, 1);

options = multigpOptions('ftc');
options.oType = {'conti', 'binary'};
options.M = M;
options.q = d;
options.params = params;
options.nlf = 1;  % no. of latent function
options.train = 0;
options.bias = bias;
    
for loop = 1:T

    rng('shuffle');
    S = 0;
    for i=1:M
        tmp = start_cnn{loop};
        initX{i} = tmp.X{i};
        initY{i} = tmp.y{i};
    end

    [model, Xobs, Yobs] = initializeBO(options, 'cnn', task, xmin, xmax, initX, initY);   
    
    [result.umax_f(loop,1), result.umax_x(loop, 1, :)] = getMaxMean(model, xmin, xmax, task, t, opts);
    result.time_real(loop, 1) = S;

    it = 2;    
    while S <= Budget
        
        disp(['%%%%%%%%%%%%%%%%%%%%% Step: ', num2str(loop), ', ', num2str(it),  ', cost: ', num2str(S),' %%%%%%%%%%%%%%%%%%%%%']);
        tic
        [xstar, values, params] = sampleMixedXstar(model, nSample, nFeatures, xmin, xmax);       
        target_xstar = reshape(xstar(1, :, :), nSample, d);

		EPc = getFactor_EP(target_xstar, values, params, model.bias(model.nlf+1:end));

        % The first EP step
        [fi_xstar, modelnew] = EPapproximationMixed(target_xstar, model, 100, EPc); % Eq. (5.16)
        disp('First EP approximation finished.');

        acq_func = cell(M, 1);
        entropy_given_xstar = cell(M, 1);
        c = cell(M, 1);
		
        for i=1:M
            c{i} = getFactor(values, fi_xstar.m(:, i), i, t, 'avg', model, xstar);
            entropy_given_xstar{i} = @(x) EPapproximationMixed2(x, i, fi_xstar.m(:, i), fi_xstar.v(:, i), modelnew, c{i});    
            acq_func{i} = @(x) (entropy_cgp(model, x', i) - entropy_given_xstar{i}(x'))/getCost_cnn(x, i);
        end
        t1 = toc;
        tic
        [optimum, type] = globalOptimizationNoGradient(acq_func, xmin', xmax', opts);
        t2 = toc;
        
		disp(['Next point selected: Elapsed time is ', num2str(t2), ' seconds']);
        
        optimum = optimum';
        optimumX = transX(optimum, 'cnn', true);
		
	    [ys, yt] = getObsValue_cnn_bos(optimumX, type, fname, gpu);
		
		tic
        [Xnew, ynew, noise] = updateMixedXY(model, optimum, ys, type);

        model = updateBO_EP(model, Xnew, ynew, noise);
        t3=toc;
        
        result.X(loop, it, :) = optimumX;
        result.y(loop, it) = ys;
        result.type(loop, it) = type;
        result.time_real(loop, it) = result.time_real(loop, it-1) + yt + (t1+t2+t3)/3600;   % CNN training time + BO running time
        
        % umax_f is the simulated function value.
        [result.umax_f(loop, it), result.umax_x(loop, it, :)] = getMaxMean1(model, xmin, xmax, task, t, opts, result.umax_x(loop, it-1, :));
        
        % The true function value can be achieved via training the CNN model with umax_x as the hyper-parameters
        % Only compute the true_func when it is necessary since it takes long time.
        % result.true_func(loop, it) = getObsValue_cnn_bos(result.umax_x(loop, it, :), 1, fname, gpu);
        
        S = result.time_real(loop, it);
        it = it + 1;
  
		save([Name, experimentNo, '_result'], 'result', 'model');
    end
    
end


