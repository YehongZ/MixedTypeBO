% Input Arguement:
% prior: prior.m -- \mu_0 in eq. 28
%        prior.var -- \Sigma_0 in eq. 28
% Output: f -- p(f^*|D_n, x^*)
function [ f ] = runEP(prior, ymax, sigma2_n, M, max_iter, c)
    
    prior.var = diag(prior.var);
    
    if size(prior.m, 1) ~= M 
        error('Dimension dismatch');
    end

    prior.varinv = 1./prior.var;  % \Sigma^{-1}_0
	fTilde = struct('m', zeros(M, 1), 'vinv', zeros(M, 1));

	% We repeat until the method converges
	convergence = false;
    iter = 0;
    
    f.v = prior.var;
    f.m = prior.m;

    b_ind = find(isnan(ymax));
    c_ind = setdiff((1:M), b_ind);

    while (~convergence && iter < max_iter) %#ok<*ALIGN>
%         iter
		fOld = f;
	
		% We eliminate the contribution of the approximate factor from the posterior approximation

		vBar = ((f.v).^-1 - fTilde.vinv).^-1;
        mBar = vBar .* (((f.v).^-1).*f.m - fTilde.vinv.*fTilde.m);

		% We refine the constraint factors
        mOldAux = mBar + c;
		vOldAux = vBar;
      
		mOldAux(c_ind) = mOldAux(c_ind) - ymax(c_ind);
		vOldAux(c_ind) = vOldAux(c_ind) + sigma2_n(c_ind);

		alpha = mOldAux ./ sqrt(vOldAux); 

		ratio = exp((-0.5 * log(2 * pi) - 0.5 * alpha.^2) - logphi(alpha));
		beta = ratio .* (alpha + ratio) ./ vOldAux;
		kappa = (alpha + ratio) ./ sqrt(vOldAux);
 
		vNewinv = beta ./ (1 - beta .* vBar);
		vNewinv( abs(vNewinv) < 1e-30 ) = 1e-30;
        mNew = mBar + 1 ./ kappa;
        
		vNewinv( vBar < 0 ) = fTilde.vinv( vBar < 0 );
		mNew( vBar < 0 ) = fTilde.m( vBar < 0 );

		failedEPupdate = true;
        damping = 1;
		while (failedEPupdate) 

			% We do damping
			mNew = mNew * damping + fTilde.m * (1 - damping);
			vNewinv = vNewinv * damping + fTilde.vinv * (1 - damping);

			% We verify that the posterior is well defined if not we increase the damping

			[ ~, eigenvalues ] = eig(diag(vNewinv +  prior.varinv));
            
			%if (any(1 ./ diag(eigenvalues) <= 1e-30))
            %    damping = damping * 0.8;
            %    disp(num2str(damping))
            %else
				failedEPupdate = false;
			%end
		end

		fTilde.m = mNew;
		fTilde.vinv = vNewinv;

		% We update the posterior marginals

		f.v = 1./(prior.varinv + fTilde.vinv);
        f.m = f.v .* (fTilde.vinv.*fTilde.m + prior.m./prior.var);
    
		change = max(abs(f.m ./ f.v - fOld.m ./ fOld.v));
		change = max(change, max(abs(f.v.^-1 - fOld.v.^-1)));

		if (change < 1e-6)
			convergence = true;
		end
	
% 		damping = damping * 0.99;

%		fprintf(1, '%f\n', change);
        iter = iter + 1;
    end

