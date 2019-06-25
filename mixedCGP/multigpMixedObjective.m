function [f, model] = multigpMixedObjective(params, model)

% Wrapper function for multigpMixedOptimise

model = modelExpandParam(model, params);
model.params = params;

if find(strcmp(model.oType, 'binary'))
    model = multigpMixedEP(model, 20);
    f = - multigpMixedLogLikelihood(model);
else
    f = - multigpLogLikelihood(model);
end
