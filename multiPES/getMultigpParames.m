% Pinv: nlf*d
% Ptinv: M*d
% sigma: nlf*M
% sigma2_n: 1*M
function [Pinv, Ptinv, sigma, sigma2_n, sigma_raw] = getMultigpParames(model)

M = model.M;
d = model.q;
nlf = model.nlf;
params = modelExtractParam(model);

Pinv = zeros(nlf, d);
for i=1:nlf
    index = paramNameRegularExpressionLookup(model,['multi ', num2str(i),' .* ', num2str(i) ...
                ' inverse width latent ']);
    Pinv(i, :) = 1./exp(params(index));
end

Ptinv = zeros(M, d);
for i=1:M
    index = paramNameRegularExpressionLookup(model,['multi .* ' num2str(i) ...
                ' inverse width output']);
    Ptinv(i, :) = 1./exp(params(index));
end

sigma = zeros(nlf, M);
sigma_raw = zeros(nlf, M);
for i=1:nlf

    index = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
                ' .* sensitivity']);

    for j=1:M
%        P = 2.*Ptinv(j, :) + Pinv(i, :);
        kern_norm = 1; % sqrt(sqrt(prod(Pinv(i, :))/prod(P)));
        sigma(i, j) = kern_norm.*params(index(j));
        sigma_raw(i, j) = params(index(j));
    end

end

sigma2_n = exp(params(end-M+1:end));
%sigma2_n(sigma2_n < 0.0001) = 0.0001;
