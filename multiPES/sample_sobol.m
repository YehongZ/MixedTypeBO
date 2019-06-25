function [ret] = sample_sobol(dim, num, xmin, xmax)
        
    p = sobolset(dim, 'Skip', 1000, 'Leap', 100);
    p = scramble(p,'MatousekAffineOwen');
    
    ret = net(p, num);
    
    scale = xmax - xmin;
    ret = ret.*repmat(scale, num, 1) + repmat(xmin, num, 1);
