% Sample x* using mixed-type random features method
function [samples, fvals, params] = sampleMixedXstar(model, nSamples, nFeatures, xmin, xmax)

    d = model.q;    % input dimension
    N = model.N;
    y = model.m;
    nlf = model.nlf;    % no. of latent function
    M = model.M;
    
    if nargout == 3
        params = cell(nSamples, 1);
    end

    samples = zeros(M, nSamples, d);
    fvals = zeros(nSamples, M);
    
    [Pinv, Ptinv, sigma, sigma2_n] = getMultigpParames(model);

    Sigma_n = [];
    for i=1:M
        switch model.oType{i}
            case 'conti'
                Sigma_n = [Sigma_n, repmat(sigma2_n(i), 1, size(model.X{i+nlf}, 1))];
            case 'binary'
                Sigma_n = [Sigma_n, model.noise{i}'];
            otherwise
                error('Invalid output type!')
        end
    end
    Sigma_n = diag(Sigma_n);
%     imagesc(Sigma_n)
    
    for i = 1:nSamples %#ok<*ALIGN>

        W_all = [];
        b_all = [];
        noise = randn(nFeatures*nlf, 1);
        trans_all = [];

        for q=1:nlf
            % draw the random features
            W = randn(nFeatures, d) ./ repmat(sqrt(Pinv(q, :)), nFeatures, 1);
            b = 2 * pi * rand(nFeatures, 1);
            
            W_all = [W_all; W];
            b_all = [b_all; b];

            trans = zeros(nFeatures, M);
            for j = 1:M         
                trans(:, j) = (sigma(q, j)*sqrt(2/nFeatures)) .* diag(exp(-0.5 .* W * diag(Ptinv(j, :)) * W'));
            end
            trans_all = [trans_all; trans];
        end

        % We sample from the posterior on the coefficients
        if (N > 0)
            Phi = [];
            for j = 1:M
                X = model.X{j+nlf};
                phi = diag(trans_all(:, j)) * cos(W_all * X' + repmat(b_all, 1, size(X, 1)));
                Phi = [Phi, phi];
            end

            % N<m, m*m matrix (A=Phi * \Sigma_n * Phi^\top + I) is rank N -> singular
            if (N < nFeatures*nlf)
                A = Sigma_n + Phi' * Phi;
                m = Phi*(A\y);
                theta = m + chol(eye(nFeatures*nlf) - Phi * (A\Phi'))' * mean(noise, 2);
            else
                A = Phi * (Sigma_n \ Phi') + eye(nFeatures*nlf);
                m = (A \ Phi) * (Sigma_n \ y);
                L = chol(A);
                theta = m + L \ mean(noise, 2);
            end
        else
            theta = mean(noise, 2);
        end
        
        for j=1:M
            % We specify the objective function and its gradient
            % x need to be a column vector
            % negative: minimization -> maximization
            targetVector = @(x) (theta' * diag(trans_all(:, j)) * cos(W_all * x + repmat(b_all, 1, size(x, 2))))' + model.bias(j+model.nlf);
            targetVectorGradient = @(x) -(theta' * diag(trans_all(:, j))) * (repmat(sin(W_all * x + b_all), 1, d) .* W_all);

            % We do global optimization
            opts.discrete = 1;
            opts.dis_num = 10;
            opts.direct.showits = 0;

            [sample, fval]= globalOptimization(targetVector, targetVectorGradient, xmin', xmax', opts);
            
            samples(j, i, :) = sample;
            fvals(i, j) = -fval;
        end
        
        if nargout == 3
            params{i}.theta = theta;
            params{i}.trans = trans_all;
            params{i}.W = W_all;
            params{i}.b = b_all;
        end
    end
    


