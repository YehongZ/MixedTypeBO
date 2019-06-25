function [K] = kernComputeFastOne(x1, type1, x2, type2, params)

    nlf = size(params.sigma, 1);
    d = size(params.Pinv, 2);
    
    K = zeros(size(x1, 1), 1);
    
    Ankinv = params.Ptinv(type1, :);
    Amkinv = params.Ptinv(type2, :);
    
    for i=1:nlf
        
        Bkinv = params.Pinv(i, :);
        sigma1_y = params.sigma(i, type1);
        sigma2_y = params.sigma(i, type2);
        
        mu_n = zeros(1, d);
        mu_m = zeros(1, d);

        P = Ankinv + Amkinv + Bkinv;
        ldet = prod(P)^(1/2);
        Linv = sqrt(1./P);
        Linvx = (x1- repmat(mu_n,size(x1,1),1))*diag(Linv);
        Linvx2 = (x2- repmat(mu_m,size(x2,1),1))*diag(Linv);
        n2 = dist2(Linvx, Linvx2);
        
        fNumPqr = prod(2*Ankinv+ Bkinv)^(1/4);
        fNumPsr = prod(2*Amkinv + Bkinv)^(1/4);
        K = K + sigma1_y*sigma2_y*fNumPqr*fNumPsr/ldet*exp(-0.5*n2);
    end



