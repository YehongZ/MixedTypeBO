function [param, names] = multigpMixedExtractParam(model)

% Extract the parameters of a multigpMixed model.
% Output:
%   param: a vector of parameters from the model.
%   names: cell array of parameter names.
% Input:
%   model: the model structure 
%
% Copyright (c) 2008 Mauricio A Alvarez
% Copyright (c) 2008 Neil D. Lawrence
% This function is modified from MULTIGPEXTRACTPARAM by Yehong for fitting the
% mixed-type CMOGP algorithm.

if nargout>1
    [kernParams, kernParamNames] = kernExtractParam(model.kern);  
else
    kernParams = kernExtractParam(model.kern);
end

kernParams = real(kernParams);

% Check if the output scales are being learnt.
if model.learnScales
    fhandle = str2func([model.scaleTransform 'Transform']);
    scaleParams = fhandle(model.scale, 'xtoa');
    if returnNames
        for i = 1:length(scaleParams)
            scaleParamNames{i} = ['Output Scale ' num2str(i)];
        end
    end
else
    scaleParams = [];
    scaleParamNames = {};
end

% Check if there is a mean function.
if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
    if nargout>1
        [meanFuncParams, meanFuncParamNames] = meanExtractParam(model.meanFunction);
        for i = 1:length(meanFuncParamNames)
            meanFuncParamNames{i} = ['mean Func ' meanFuncParamNames{i}];
        end
    else
        meanFuncParams = meanExtractParam(model.meanFunction);
    end
else
    meanFuncParamNames = {};
    meanFuncParams =[];
end

% Check if there is a parameter beta

if isfield(model, 'beta') && ~isempty(model.beta)
    fhandle = str2func([model.betaTransform 'Transform']);
    betaParams = fhandle(model.beta, 'xtoa');    
    if nargout>1
        betaParamNames = cell(model.nout,1);
        for i = 1:length(betaParams)
            betaParamNames{i} = ['Beta ' num2str(i)];
        end        
    end
else
    betaParamNames = {};
    betaParams =[];
end

% Check if there is a noise model
paramPart = [kernParams scaleParams meanFuncParams betaParams];
if nargout > 1
    names = {kernParamNames{:}, meanFuncParamNames{:}, ...
        scaleParamNames{:}, betaParamNames{:}};
end

if isfield(model, 'fix')
    for i = 1:length(model.fix)
        paramPart(model.fix(i).index) = model.fix(i).value;
    end
end


switch model.approx
    case 'ftc'
        param = paramPart;
    case {'dtc','fitc','pitc'}
        if nargout>1
            [param, names] = spmultigpExtractParam(model, paramPart, names);
        else
            param = spmultigpExtractParam(model, paramPart);
        end
    otherwise
        error('Unknown approximation')
end
