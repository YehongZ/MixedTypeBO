function  model = multigpMixedCreate(q, d, X, y, options)
% Creates a mixed-type CMOGP model.
% Returns:
%   MODEL: the structure for the multigpMixed model
% Arguments:
%   q: input dimension.
%   d: output dimension size.
%	y: a cell of observed outputs
%	X: a set of inputs
%	options: contains the options for the multigpMixed model
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
% MODIFICATIONS : David Luengo, 2009
% This function is modified from MULTIGPCREATE by Yehong for fitting the
% mixed-type CMOGP algorithm.	

if iscell(X)
    if size(X{end}, 2) ~= q
        error(['Input matrix X does not have dimension ' num2str(q)]);
    end
else
    if size(X, 2) ~= q
        error(['Input matrix X does not have dimension ' num2str(q)]);
    end
end

if iscell(y)
    % ensure it is a row vector of cells.
    y = y(:)';
    if size(y, 2) ~= d
        error(['Target cell array Y does not have dimension ' num2str(d)]);
    end
    for i = 1:size(y, 2)
        if(size(y{i}, 2)>1)
            error('Each element of the cell array should be a column vector.')
        end
    end
else
    if size(y, 2)~=d
        error(['Target matrix Y does not have dimension ' num2str(d)]);
    end
end


model.type = 'multigpMixed';
model.q = q;
model.M = options.M;

% Short hand for the kernel type
model.kernType = options.kernType;

% Number of latent functions
model.nlf = options.nlf;
if isfield(options, 'typeLf') && ~isempty(options.typeLf)
   model.typeLf = options.typeLf; 
end

if isfield(options, 'oType') && ~isempty(options.oType)
   model.oType = options.oType; 
end

% Number of output functions
model.d = d;

switch options.approx
    case 'ftc'
        model.nout = d - options.nlf;
    case {'dtc','fitc','pitc'}
        error('Not implement yet.')
end

model.approx = options.approx;

% Set up default scale and bias for outputs
if isfield(options, 'scale') && ~isempty(options.scale)
    model.scale = options.scale;
else
    model.scale = ones(1, model.d);
end
if isfield(options, 'bias') && ~isempty(options.bias)
    model.bias = options.bias;
else
    model.bias = zeros(1, model.d);
end

% Initialization of the model
switch model.approx
    case 'ftc'
        model.X = X;
        model.oY = y;
        model.y = [];
        
        for i = 1:length(y)-model.nlf
            switch model.oType{i}
                case 'conti'
                    model.y = [model.y; y{i+model.nlf}];
                case 'binary'
                    model.y = [model.y; y{i+model.nlf}];%zeros(length(y{i+model.nlf}), 1)];
                otherwise
                    error('Invalid output type!');
            end
        end
        model.N = size(model.y,1);
    case {'dtc','fitc','pitc'}
        
    otherwise
        error('Unknown model approximation')
end

model.includeInd = options.includeInd;	
model.tieIndices = options.tieOptions.tieIndices;   
model.includeNoise = options.includeNoise;  % 1-include noise w(x)
kernType = cell(1,sum(options.nlf) + options.includeNoise + ...
    options.includeInd);
cont = 0;

% This checks if there are different clases of latente forces
if isfield(options, 'typeLf') && ~isempty(options.typeLf)
    if options.nlf ~= sum(options.typeLf) 
        error('The number of latent functions per type is different to the total number of latent functions')
    end
    for i = 1:options.typeLf(1)
        cont = cont + 1;
        kernType{cont} = multigpKernComposer(options.kernType, model.d, model.nlf, model.approx, i);
    end    
    for i = 1:options.typeLf(2)
        cont = cont + 1;
        kernType{cont} = multigpKernComposer([options.kernType 'white'], model.d, model.nlf, model.approx, i+options.typeLf(1));
    end    
else
    for i = 1:options.nlf
        cont = cont + 1;
        kernType{cont} = multigpKernComposer(options.kernType, model.d, model.nlf, model.approx, i, options);
    end
end
% To include independent kernel
if model.includeInd
    cont = cont + 1;
    kernType{cont} = multigpKernComposer('rbf', model.d, model.nlf, model.approx);
end
% To include noise
if model.includeNoise
    cont = cont + 1;
    kernType{cont} = multigpKernComposer('white', model.d, model.nlf, model.approx);
end
model.kern = kernCreate(X, {'cmpnd', kernType{:}}); %cmpnd - sum of the individual kernels (f_q(x)+w_q(x))

if ~isfield(options, 'isSpeedUp') || (isfield(options, 'isSpeedUp') && options.isSpeedUp ~= 2)    
    % These are additional options for some kernels that must be handled outside
    % kernCreate. For example, for the GG kernel, the user can choose ARD or
    % not.
    fhandle = [model.kernType 'MultigpKernOptions'];
    if exist(fhandle, 'file')
        fhandle = str2func(fhandle);
        model = fhandle(model, options);
    end
end

if isfield(options, 'optimiser') && ~isempty(options.optimiser)
    model.optimiser = options.optimiser;    
end

% NEIL: Learning of scales hasn't been included, although it should be.
model.learnScales = options.learnScales;
model.scaleTransform = optimiDefaultConstraint('positive');

model.nParams = 0;

% Extract top level parameters from kernels.
model = multigpUpdateTopLevelParams(model); %%useless in our multi-output gp models

% Set up a mean function if one is given.
if isfield(options, 'meanFunction') && ~isempty(options.meanFunction)
    if isstruct(options.meanFunction)
        model.meanFunction = options.meanFunction;
    else
        if ~isempty(options.meanFunction)
            model.meanFunction = meanCreate(q, model.nout, X, y, options.meanFunctionOptions);
        end
    end
    model.nParams = model.nParams + model.meanFunction.nParams;
end

% Tie options according to the particular kernel employed

%%% Yehong: tieInd - groups of parameters of a model to be seen as one parameter during optimisation of the model 
%%% Yehong: e.g. sigma2_u/precision_u in gaussian and gg kernel

fhandle = [model.kernType 'MultigpTieParam'];
if exist(fhandle, 'file')
    fhandle = str2func(fhandle);
    tieInd = fhandle(model, options);
else
    if ~isfield(options, 'isSpeedUp') || (isfield(options, 'isSpeedUp') && options.isSpeedUp ~= 2)
        error('Function for tying parameters for this kernel type not implemented yet')
    else
        warning('multigpCreate:NoTyingFunctionGlobalKernel','Function for tying parameters for this kernel type not implemented yet')
    end
end

model.nParams = model.nParams + model.kern.nParams;

% Fix parameters options according to the particular kernel employed
fhandle = [model.kernType 'MultigpMixedFixParam'];
if exist(fhandle, 'file')
    fhandle = str2func(fhandle);
    model = fhandle(model, options);
% else
%     warning('multigp:FixedParam','Function for fixing parameters for this kernel type not implemented yet')
end
if ~isfield(options, 'isSpeedUp') || (isfield(options, 'isSpeedUp') && options.isSpeedUp ~= 2)
    if isfield(options, 'fixInducing') && ~options.fixInducing && ...
            isfield(options, 'tieInducing') && ~isempty(options.tieInducing) && ...
            options.tieInducing && options.nlf > 1
        tieInd(end+1:end+(model.k(1)*model.q)) = spmultigpTiePseudoInputs(model);
    end
    model = modelTieParam(model, tieInd);
else
    if exist('tieInd', 'var') && ~isempty(tieInd)
        model = spmultimodelTieParam(model, tieInd);
    end
end

if isfield(options, 'useKernDiagGradient') && options.useKernDiagGradient
    model.useKernDiagGradient = options.useKernDiagGradient;    
end

%model.alpha = [];

model.noise = cell(model.M, 1);
if strcmp(model.kern.type, 'cmpnd') == 1 && length(model.kern.comp) == 2
    for i=1:model.M
        switch model.oType{i}
            case 'conti'
                model.noise{i} = model.kern.comp{2}.comp{i+model.nlf}.variance;
            case 'binary'
                %model.kern.comp{2}.comp{i+model.nlf}.variance = 0;
                model.noise{i} = exp(-2)*ones(size(model.X{i+model.nlf}, 1), 1);
            otherwise
                error('Invalid output type!');
        end
    end
else
    error('Invalid kernel type!');
end

params = modelExtractParam(model);
model = modelExpandParam(model, params);
model.params = params;

end
