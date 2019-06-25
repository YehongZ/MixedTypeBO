function [gParam, gX_u] = multigpMixedLogLikeGradients(model)

% Compute the gradients for the parameters and X_u.
%
% COPYRIGHT : Mauricio A Alvarez, 2008, 2010
% This function is modified from MULTIGPLOGLIKEGRADIENTS by Yehong for fitting the
% mixed-type CMOGP algorithm.

gX_u = [];
g_scaleBias = gpScaleBiasGradient(model);

switch model.approx
    case 'ftc'
        covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;
        covGrad = 0.5*covGrad;        
        if sum(cellfun('length',model.X)) == model.N,
            base = 0;
        else
            base = model.nlf;
        end 
        covGrad2 = zeros(base + model.N);
        covGrad2(base+1:base+model.N,base+1:base+model.N) = covGrad;
        covGrad = covGrad2;
        gParam = kernGradient(model.kern, model.X, covGrad);
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
             if  strcmp(model.kernType,'lfm') || strcmp(model.kernType,'sim') ...
                    || strcmp(model.kernType,'lfmwhite') ...
                    || strcmp(model.kernType,'simwhite') ...
                    || strcmp(model.kernType,'simnorm') 
                gmuFull = model.m'*model.invK;
                gmu = zeros(1,model.d);
                startVal = 1;
                endVal = 0;
                for j=1:model.d-base,
                    endVal =  endVal + length(model.X{j+base});
                    gmu(j+base) = sum(gmuFull(startVal:endVal));
                    startVal = endVal + 1;
                end
                gmu = gmu(model.nlf+1:end);
                g_meanFunc = meanGradient(model.meanFunction, gmu);
            else
                g_meanFunc = [];
            end
        else
            g_meanFunc = [];
        end
        gParam = [gParam g_scaleBias g_meanFunc];
        if isfield(model, 'fix')
            for i = 1:length(model.fix)
                gParam(model.fix(i).index) = 0;
            end
        end
    case {'dtc', 'fitc', 'pitc'}
        % Sparse approximations.
        if isfield(model, 'beta') && ~isempty(model.beta)
            [dLdKyy, dLdKuy, dLdKuu, dLdmu, gBeta] = spmultigpLocalCovGradient(model);
            fhandle = str2func([model.betaTransform 'Transform']);
            gBeta = gBeta.*fhandle(model.beta, 'gradfact');
        else
            [dLdKyy, dLdKuy, dLdKuu, dLdmu] = spmultigpLocalCovGradient(model);
            gBeta = [];
        end
        if ~model.fixInducing
            [gParam, gX_u] = sparseKernGradient(model, dLdKyy, dLdKuy, dLdKuu);
        else
            gParam = sparseKernGradient(model, dLdKyy, dLdKuy, dLdKuu);
        end
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
            if  strcmp(model.kernType,'lfm') || strcmp(model.kernType,'sim') ...
                    || strcmp(model.kernType,'lfmwhite') ...
                    || strcmp(model.kernType,'simwhite') ...
                    || strcmp(model.kernType,'simnorm')
                gmu = zeros(1,model.nout);
                for j=1:model.nout,
                    gmu(j) = sum(dLdmu{j});
                end
                g_meanFunc = meanGradient(model.meanFunction, gmu);
            else
                g_meanFunc = [];
            end
        else
            g_meanFunc = [];
        end
        gParam = [gParam g_scaleBias g_meanFunc gBeta];
        if isfield(model, 'fix')
            for i = 1:length(model.fix)
                gParam(model.fix(i).index) = 0;
            end
        end
        % if there is only one output argument, pack gX_u and gParam into it.
        if ~model.fixInducing || nargout > 1
            gParam = [gParam gX_u(:)'];
        end
end





