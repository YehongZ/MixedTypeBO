function [ret] = sample_rand(dim, num, xmin, xmax)
        
    ret = rand(num, dim);
    
    scale = xmax - xmin;
    ret = ret.*repmat(scale, num, 1) + repmat(xmin, num, 1);
