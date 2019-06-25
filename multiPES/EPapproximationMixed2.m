% EP approximation: Approximating terms that depend on x
function [ret] = EPapproximationMixed2(x, type, fi_xstar_m, fi_xstar_v, modelnew, c)
    
    nSamples = size(fi_xstar_m, 1);
    
    model = modelnew{1,1};
    sigma2_n = getNoise(model);   
    a = [-1,1];
    
    nlf = model.nlf;
    s = 1+nlf;
    t = type+nlf;
    
    ret = 0;
    no = 0;
    for i=1:nSamples
        
        fplus = fplusComputeOne(x, type, fi_xstar_m(i), fi_xstar_v(i), modelnew{type,i}, s, t);
        m = a*fplus.m;
        v = a*fplus.v*a';

        if v < 0
%             warning('Negative v. Your posterior will equal to the prior.');
            vinv = 0;
        else
            vinv = 1/v;
        end
        
        alpha = (c-m)*sqrt(vinv);
        gamma = exp((-0.5 * log(2 * pi) - 0.5 * alpha.^2) - logphi(alpha));

        varNew = fplus.v - vinv * gamma * (gamma-(m-c)*sqrt(vinv))* fplus.v * (a'*a) * fplus.v;
        muNew = fplus.m - gamma*sqrt(vinv)*(fplus.v*a');        
        res = min(varNew(2,2), fplus.v(2,2));
        if res < varNew(2,2)
            muNew = fplus.m;
        end
        switch model.oType{type}
            case 'conti'
                res = res + sigma2_n(type);
                if res <= 0
                    no = no+1;
                else
                    ret = ret + 0.5*log(res);
                end
            case 'binary'
                if 1+res <= 0
                    no = no+1;
                else
                    p = exp(logphi(muNew(2)/sqrt(1+res)));
                    h = - p*log(p) - (1-p)*log(1-p);
                    if isnan(h)
                        h = 0;
                    end
					ret = ret + h;
                end
            otherwise
                error('Invalid output type!');
        end
    end
    
    if no > 0
    	disp(['The predictive variance given x* is negative, ', num2str(no), ' results ignored']);
    end
    ret = ret/(nSamples-no);
end

function fplus = fplusComputeOne(x, type, fi_xstar_m, fi_xstar_v, model, s, t)

    M = model.M;
    nlf = model.nlf;
    [mu, v, KX_star] = multigpPosteriorSinglef_fast(model, x, type);

    nObs = zeros(1, M+nlf);
    for j=s:s+M-1
        nObs(j) = size(model.X{j}, 1);
    end
    
    Kinvk = KX_star'*model.invK;
    m1 = Kinvk(sum(nObs(1:t)));
    
    fplus.m = [fi_xstar_m; mu];
    fplus.v = [fi_xstar_v, fi_xstar_v*m1; fi_xstar_v*m1, v+m1^2*fi_xstar_v]; 

end

