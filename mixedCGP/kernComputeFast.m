function [ret] = kernComputeFast(x1, x2, type, params)
    
    [nlf, M] = size(params.sigma);
    
    if iscell(x1)
        N = 0;
        for i=1:M
            N = N + size(x1{nlf+i}, 1);
        end

        ret = [];
        for i=1:M
            x = x1{nlf+i};
            ret = [ret; kernComputeFastOne(x, i, x2, type, params)];
        end
    else
        ret = kernComputeFastOne(x1, type, x2, type, params);
    end
