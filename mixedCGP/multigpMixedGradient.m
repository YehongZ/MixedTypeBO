function g = multigpMixedGradient(params, model, update)

% Gradient wrapper for a multigpMixed model.

model = modelExpandParam(model, params);
model.params = params;

if find(strcmp(model.oType, 'binary'))
    if nargin == 3 && update
        model = multigpMixedEP(model, 20);
    end
end

g = - modelLogLikeGradients(model);