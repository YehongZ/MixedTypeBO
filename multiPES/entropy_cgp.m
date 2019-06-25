% Compute the predictive entropy given a multigpMixed model  
function [ret] = entropy_cgp(model, x, type)  

    [~, var] = multigpMixedPosteriorMeanVar(model, x);
    
    switch model.oType{type}
        case 'conti'
            ret = 0.5*log(var{type+model.nlf});
        case 'binary'
            p = var{type+model.nlf};
            ret = -p.*log(p)-(1-p).*log(1-p);
            ret(isnan(ret)) = 0;
        otherwise
            error('Invalid output type!');
    end
    