function [mu, varSig, ret_kx_star] = multigpPosteriorSinglef(model, X, type)

% Compute the mean and variance of the posterior distribution without noise.
% Output:
%	  mu: vector containing the posterior mean
%	  varSig: posterior variance matrix
% Input:
%	  model: the model for which posterior is to be computed.
%	  X: an input to be predicted
%	  type: type of ouptut (default 0: predict for all types of output)

if nargin<3
    type = 0;
end

% If X is a vector assume it applies for all outputs.
if ~iscell(X)
    nX = size(X, 1);
    xtemp = X;
    X = cell(1, model.d);
    for i = 1:model.d
        X{i} = xtemp;
    end
    
    % Yehong: pitc
    X1 = cell(1, model.d+model.nlf);
    for i = 1:model.d+model.nlf
        X1{i} = xtemp;
    end
end

model.kern.type = 'cmpndno';

switch model.approx
    case 'ftc'
        KX_star = kernCompute(model.kern, model.X, X); %(M+nlf)*(M*|X|)
        Kstar = kernCompute(model.kern, X);     %((M+nlf)*|X|)*((M+nlf)*|X|)
 
        % Again, this is a bit of a hack to allow the inclusion of the
        % latent structure in the kernel structure.
        if sum(cellfun('length',model.X)) == model.N,
            base = 0;
        else
            base = model.nlf;
        end
        
        if type == 0
            KX_star = KX_star(base+1:end,:);    % N*(M*|X|)
            muTemp = KX_star'*model.alpha;      % (M+1)*1
            %diagK = kernDiagCompute(model.kern, X);
            varsigTemp = Kstar - KX_star' * model.invK * KX_star;

            muTemp = (muTemp' .* model.scale + model.bias)';
            varsigTemp = diag(model.scale)*varsigTemp*diag(model.scale);

            mu = muTemp(base+nX:end);
            varSig = varsigTemp(base+nX:end, base+nX:end);

            ret_kx_star = KX_star;
            
        else

            KX_star = KX_star(base+1:end, nX*(base+type-1)+1:nX*(base+type-1)+nX);    % N*|X|
            muTemp = KX_star'*model.alpha;      % |X|*1
            
            Kstar = Kstar(nX*(base+type-1)+1:nX*(base+type-1)+nX, nX*(base+type-1)+1:nX*(base+type-1)+nX); % |X|*|X|
            varsigTemp = Kstar - KX_star' * model.invK * KX_star;

            mu = muTemp * model.scale(type+base) + model.bias(type+base);
            varSig = model.scale(type+base).^2 * varsigTemp;

            ret_kx_star = KX_star;
        end
        
    otherwise
        error('multigpPosteriorSinglef not implemented yet');
end





