% First EP approximation step (independent of x) 
function [ret, modelnew] = EPapproximationMixed(xstar, model, max_iter, c)
    
    nSamples = size(xstar, 1);
	f_star = cell(nSamples, 1);
    
    % p(f_*|D_n)
    for i = 1 : nSamples
        [f_star{i}.m, f_star{i}.var] = multigpPosteriorSinglef(model, xstar(i, :));
    end
    
    M = model.M;
    nlf = model.nlf;

    ymax = zeros(M, 1);
    sigma2_n = getNoise(model);
    
    for i=1:M
        if strcmp(model.oType{i}, 'conti') == 1
            y = model.oY{i+nlf};
            ymax(i) = max(y);
        else
            ymax(i) = nan;
        end
    end
    
    M = model.M;
    ret.m = zeros(nSamples, M);
    ret.v = zeros(nSamples, M);
    for i = 1 : nSamples
        tmp_ret = runEP(f_star{i}, ymax, sigma2_n, M, max_iter, c);
        ret.m(i, :) = tmp_ret.m;
        ret.v(i, :) = tmp_ret.v;    
    end

    modelnew = cell(M, nSamples);
    options = model.options;
    
    % Construct the new models given (x*, f*_i)
    for i=1:M
        for j=1:nSamples        
            [Xnew, ynew] = updateXY(model, xstar(j, :), ret.m(j, i), i);
            tmodel = model;
            if strcmp(model.oType{i}, 'binary') == 1    
                tmodel.noise{i} = [tmodel.noise{i}; 0];
            end
            modelnew{i, j} = multigpMixedUpdate(tmodel, Xnew, ynew, options);
            
            [Pinv, Ptinv, ~, sigma2_n, sigma] = getMultigpParames(modelnew{i, j});
            modelnew{i, j}.hypers.Pinv = Pinv;
            modelnew{i, j}.hypers.Ptinv = Ptinv;
            modelnew{i, j}.hypers.sigma = sigma;
            modelnew{i, j}.hypers.sigma2_n = sigma2_n;
        end
    end
    
    
    
    
    
    
    
