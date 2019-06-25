% Extract the parameters of the binary output
% Outputs:
% K: K_{X_B, X_B}
function [K, m, y, index, ndata] = multigpMixedExtractBinary(model)

index = [];
ndata = [];
y = [];
nX = [];
N = 0;
m = [];     
for i = 1:model.M
    if strcmp(model.oType{i}, 'binary') == 1
        index = [index, i];
        ndata = [ndata, N];
        y = [y; model.oY{i+model.nlf}];
        nX = [nX, N+1:N+length(model.oY{i+model.nlf})];
        m = [m; model.bias(i+model.nlf)*ones(length(model.oY{i+model.nlf}), 1)];
    end
    N = N + length(model.oY{i+model.nlf});
end

K = model.fK(nX, nX);
