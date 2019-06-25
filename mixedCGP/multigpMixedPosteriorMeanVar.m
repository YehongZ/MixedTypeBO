function [mu, varsig] = multigpMixedPosteriorMeanVar(model, X, computeAll)

% Compute the posterior mean and variance given a multigpMixed model
% Input:
%   model : the model for which posterior is to be computed.
%	X : cell array containing locations where outputs are to be computed.
% Output:
%	mu : cell array containing mean posterior vectors.
%	varsig : cell array containing the variance posterior vector
%
% COPYRIGHT : Mauricio A Alvarez, 2008, 2009, 2010
% COPYRIGHT : Neil D. Lawrence, 2008
% MODIFICATIONS : David Luengo, 2009
% MODIFICATIONS : Mauricio A. Alvarez, 2010
% This function is modified from MULTIGPPOSTERIORMEANVAR by Yehong for fitting the
% mixed-type CMOGP algorithm.

if nargin<3
    computeAll = true;
end

% If X is a vector assume it applies for all outputs.
if ~iscell(X)
    xtemp = X;
    X = cell(1, model.d);
    for i = 1:model.d
        X{i} = xtemp;
    end
end

switch model.approx
    
    case 'ftc'
        mu = cell(model.d,1);
        varsig = cell(model.d,1);
        KX_star = kernCompute(model.kern, model.X, X);
        
        if sum(cellfun('length',model.X)) == model.N
            base = 0;
        else
            base = model.nlf;
        end
        
        KX_star = KX_star(base+1:end,:);
        muTemp = KX_star'*model.alpha;
        diagK = kernDiagCompute(model.kern, X);
        Kinvk = model.invK*KX_star;
        varsigTemp = diagK - sum(KX_star.*Kinvk, 1)';
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
            m = meanCompute(model.meanFunction, X, model.nlf);
        end
        
        startVal=1;
        endVal=0;
        for i=1:length(X)
            endVal = endVal + size(X{i},1);
            if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                mu{i}  = muTemp(startVal:endVal, 1)*model.scale(i) ...
                    + m(startVal:endVal, 1)+ model.bias(i);
            else
                mu{i}  = muTemp(startVal:endVal, 1)*model.scale(i) + model.bias(i);
            end
            varsig{i} = varsigTemp(startVal:endVal, 1)*model.scale(i)*model.scale(i);
            startVal = endVal+1;
        end

        if isfield(model, 'oType')
            for i = 1:model.M
                if strcmp(model.oType{i}, 'binary') == 1
                    p_positive = exp(logphi(mu{i+base}./sqrt(1+varsig{i+base})));
                    varsig{i+base} = p_positive;
                    mu{i+base} = sign(p_positive-0.5);
                end
            end
        end

    otherwise
        error('multigpMixedPosteriorMeanVar not yet implemented');
end


