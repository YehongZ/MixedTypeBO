
function [sigma2_n] = getNoise(model)

M = model.M;
params = model.params;

sigma2_n = exp(params(end-M+1:end)');
