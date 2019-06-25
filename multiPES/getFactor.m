% Compute the slack variables c_i
function [ret] = getFactor(fi_max, fi_xstar_mu, type, t, method, model, xstar)

    if type == t
        ret = 0;
    else
        switch method
            case 'max'
                ret = max(0, max(fi_max(:, type) - fi_xstar_mu));
            case 'avg'
                ret = mean(max(fi_max(:, type) - fi_xstar_mu, 0));
            case 'avg2'
                x = reshape(xstar(type, :, :), size(xstar, 2), size(xstar, 3));
                [mu, ~] = multigpPosteriorMeanVar(model, x);              
                ret = mean(abs(mu{type+model.nlf} - fi_xstar_mu));
            case 'avg1'
                x = reshape(xstar(type, :, :), size(xstar, 2), size(xstar, 3));
                [mu, ~] = multigpPosteriorMeanVar(model, x);              
                ret = mean(max(0, mu{type+model.nlf} - fi_xstar_mu));
            case 'fixed'
                ret = 0.2;
            case 'max_mean'
                [mu, ~] = multigpPosteriorMeanVar(model, xstar);
                m = mu{type+model.nlf};
                ret = max(0, m(type)-m(t));
            otherwise
                error('Invalid method for computing c.')
        end
    end
    
    
end
