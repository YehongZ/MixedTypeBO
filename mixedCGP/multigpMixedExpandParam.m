function model = multigpMixedExpandParam(model, params)

% Expand the given parameters into a multigpMixed structure.
% Input:
%   model : the model structure
%	params : a vector of new parameters for the model.
% Output:
%   model : the model structure with the new parameter vector.
%
% COPYRIGHT : Mauricio A. Alvarez, 2008
% COPYRIGHT : Neil D. Lawrence, 2008
% This function is modified from MULTIGPEXPANDPARAM by Yehong for fitting the
% mixed-type CMOGP algorithm.	

paramPart = real(params);

if isfield(model, 'fix')
    for i = 1:length(model.fix)
       paramPart(model.fix(i).index) = model.fix(i).value;
    end
end

startVal = 1;
endVal = model.kern.nParams;
kernParams = paramPart(startVal:endVal);

if length(kernParams) ~= model.kern.nParams
    error('Kernel Parameter vector is incorrect length');
end

model.kern = kernExpandParam(model.kern, kernParams);

% Check if there is a mean function.`
if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
    startVal = endVal + 1;
    endVal = endVal + model.meanFunction.nParams;
    model.meanFunction = meanExpandParam(model.meanFunction, ...
        paramPart(startVal:endVal));
end

% Check if there is a beta parameter.
if isfield(model, 'beta') && ~isempty(model.beta)
    startVal = endVal + 1;
    endVal = endVal + model.nout;
    fhandle = str2func([model.betaTransform 'Transform']);
    model.beta = fhandle(paramPart(startVal:endVal), 'atox');
end

switch model.approx
    case {'dtc','fitc','pitc'}
        error('Approximation method not implenent yet')
    otherwise
        %
end

% Keep the values of parameters related with the mean at the top level.
model = multigpUpdateTopLevelParams(model);
model.m = multigpComputeM(model);
model = multigpMixedUpdateKernels(model);

% Update the vector 'alpha' for computing posterior mean.
if isfield(model, 'alpha')
    model.alpha = multigpComputeAlpha(model);
end

if isfield(model, 'params')
    [Pinv, Ptinv, ~, sigma2_n, sigma] = getMultigpParames(model);
    model.hypers.Pinv = Pinv;
    model.hypers.Ptinv = Ptinv;
    model.hypers.sigma = sigma;
    model.hypers.sigma2_n = sigma2_n;
end

for i=1:model.M
    if strcmp(model.oType{i}, 'conti') == 1
        model.noise{i} = model.kern.comp{2}.comp{i+model.nlf}.variance;
    end
end
