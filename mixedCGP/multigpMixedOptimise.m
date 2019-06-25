function [model, fmin, params] = multigpMixedOptimise(model, display, iters)

% Optimise the hypers (params) for multigpMixed model.
% Output:
% 	model : the optimised model.
% 	params : the optimised parameter vector.
% Input:  
%   model : the model to be optimised.
%   display : whether or not to display while optimisation proceeds,
%	   set to 2 for the most verbose and 0 for the least verbose.
%	iters : number of iterations for the optimisation.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
% MODIFIED : Mauricio A. Alvarez, 2008
% This function is modified from MULTIGPOPTIMISE by Yehong for fitting the
% mixed-type CMOGP algorithm.

if nargin < 3
  iters = 1000;
  if nargin < 2
    display = 1;
  end
end

params = modelExtractParam(model);

options = optOptions;
if display
  options(1) = 1;
  if length(params) <= 100 && display > 1
    options(9) = 1;
  end
end
options(14) = iters;

if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('conjgrad');
end

if strcmp(func2str(optim), 'optimiMinimize')
  % Carl Rasmussen's minimize function 
  params = optim('multigpObjectiveGradient', params, options, model);
else
  % NETLAB style optimization.
  [params, options] = optim('multigpMixedObjective', params,  options, ...
                 'multigpMixedGradient', model);
end

%model = multigpExpandParam(model, params);

model = modelExpandParam(model, params);
model.params = params;
fmin = options(8);

if find(strcmp(model.oType, 'binary'))
    model = multigpMixedEP(model, 20);
end



