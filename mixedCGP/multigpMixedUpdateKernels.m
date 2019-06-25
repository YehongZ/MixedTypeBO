function model = multigpMixedUpdateKernels(model)

% Updates the kernel representations in the multigpMixed structure
%
% COPYRIGHT : Mauricio A. Alvarez, 2010
% This function is modified from MULTIGPUPDATEKERNELS by Yehong for fitting the
% mixed-type CMOGP algorithm.

switch model.approx
    case 'ftc'
        
        fK = real(kernCompute(model.kern.comp{1}, model.X));
        K = fK + real(kernCompute(model.kern.comp{2}, model.X));
        % This is a hack to allow the inclusion of the latent structure
        % inside the kernel structure
        if sum(cellfun('length',model.X)) == model.N
            base = 0;
        else
            base = model.nlf;
        end       
        
        model.fK = fK(base+1:end, base+1:end);   
        
        K = K(base+1:end, base+1:end);    
        if isfield(model, 'oType')
            ct = 0;
            for i = 1:model.M
                N = size(model.X{i+model.nlf}, 1);
                if strcmp(model.oType{i}, 'binary') == 1
                    if N ~= length(model.noise{i})
                        error('dimension not match')
                    end
                    K(ct+1:ct+N, ct+1:ct+N) = K(ct+1:ct+N, ct+1:ct+N) + diag(model.noise{i});
                end
                ct = ct + N;
            end           
        end
        model.K = K;
        
        [model.invK, U, jitter] = pdinv(K);
%         if any(jitter>1e-4)
%             fprintf('Warning: UpdateKernels added jitter of %2.4f\n', jitter)
%         end
        model.alpha = multigpComputeAlpha(model);
        model.logDetK = logdet(model.K, U);
    case {'dtc','fitc','pitc'}
        error('Invalid approx type!');
end


