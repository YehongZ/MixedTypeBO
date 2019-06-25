function [model, X, Y] = initializeBO(options, method, task, xmin, xmax, givenX, givenY)
    
    M = options.M;
    nlf = options.nlf;
    q = options.q;
    X = cell(M, 1);
    Y = cell(M, 1);
%     noise = task{M+1};
    
    switch method       
            
        case 'givenX'
            data_size = ones(1, M);
            scale = xmax-xmin;
            for i=1:M
                X{i} = givenX{i};%.*repmat(scale, data_size(i), 1) + repmat(xmin, data_size(i), 1);
            end
            
        case 'lr'
            for i=1:M
                X{i} = givenX{i};
                Y{i} = givenY{i};
            end
            
        case 'cnn'
            for i=1:M
                X{i} = givenX{i};
                Y{i} = givenY{i};
            end
            
        case 'rl_car'
            for i=1:M
                X{i} = givenX{i};
                Y{i} = givenY{i};
            end
            
        otherwise
            error('Invalid initialization method.')
    end
    
    if ~strcmp(method, 'lr') && ~strcmp(method, 'cnn') && ~strcmp(method, 'rl_car')
        for i=1:M
            Y{i} = getObsValue(task, X{i}, M, i);
        end
    end
    
    switch options.approx
        case 'ftc'
            Xtrain = cell(M+nlf, 1);
            ytrain = cell(M+nlf, 1);
            for j=1:nlf
               ytrain{j} = [];
               Xtrain{j} = zeros(1, q);  
            end
            s = nlf;
        otherwise
            error('Invalid approx method');
    end
    
    if ~isfield(options, 'ymean')
        options.ymean = zeros(M, 1);
    end
    
    for j=1:M
        Xtrain{j+s} = transX(X{j}, method);
        ytrain{j+s} = transY(Y{j}, options.ymean(j), method);
    end

    if options.train
        error('Not implement for training model!')
    else
        options.kernType = 'gg';
        options.optimiser = 'scgMixed';
        options.d = options.M + options.nlf;
        options.isArd = 1;
        
        model = multigpMixedCreate(options.q, options.d, Xtrain, ytrain, options);
        params = options.params;
        model = modelExpandParam(model, params);
        params = modelExtractParam(model);
        model.params = params;
        model.options = options;
        model.train = options.train;
    end
    
    