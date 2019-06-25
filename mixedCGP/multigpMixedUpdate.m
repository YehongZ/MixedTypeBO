function  model = multigpMixedUpdate(model, X, y, options, noise)
% Update the observations of a mixed-type CMOGP model.
% Returns:
%   model: the updated multigpMixed model
% Arguments:
%   model: old model
%	y: a cell of observed outputs
%	X: a set of inputs
%	options: contains the options for the multigpMixed model
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
% MODIFICATIONS : David Luengo, 2009
% This function is modified from MULTIGPCREATE by Yehong for fitting the
% mixed-type CMOGP algorithm.

model.type = 'multigpMixed';
model.approx = options.approx;

% Initialization of the model
switch model.approx
    case 'ftc'
        model.X = X;
        model.oY = y;
        model.y = [];

        for i = 1:length(y)-model.nlf
            model.y = [model.y; y{i+model.nlf}];
        end
        model.N = size(model.y,1);
    case {'dtc','fitc','pitc'}
        error('Not implement yet!')
    otherwise
        error('Unknown model approximation')
end

if ~isfield(model, 'noise')
    model.noise = cell(model.M, 1);
    for i=1:model.M
        switch model.oType{i}
            case 'conti'
                model.noise{i} = model.kern.comp{2}.comp{i+model.nlf}.variance;
            case 'binary'
                model.noise{i} = exp(-2)*ones(size(model.X{i+model.nlf}, 1), 1);
            otherwise
                error('Invalid output type!');
        end
    end  
end

if nargin == 5
    model.noise = noise;
end

params = modelExtractParam(model);
model = modelExpandParam(model, params);
model.params = params;

