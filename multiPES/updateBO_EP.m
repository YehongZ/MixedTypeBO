function  ret = updateBO_EP(model, X, Y, noise)
    
    if model.train
        error('Not implement yet');
    else
        options = model.options;
        if nargin == 4
            model = multigpMixedUpdate(model, X, Y, options, noise);
        else
            model = multigpMixedUpdate(model, X, Y, options);
        end
        ret = multigpMixedEP(model, 20);
%         ret = updateVarianceBinary(model);
    end
    
end