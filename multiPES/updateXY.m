function [Xnew, ynew] = updateXY(model, x, y, type)

    nlf = model.nlf;
    q = model.q;
    M = model.M;
    
    switch model.approx
        case 'ftc'
            Xnew = cell(M+nlf, 1);
            ynew = cell(M+nlf, 1);
            for j=1:nlf
               ynew{j} = [];
               Xnew{j} = zeros(1, q);  
            end
            s = 1+nlf;
            t = type+nlf;
        case {'dtc','fitc','pitc'}
            Xnew = cell(M, 1);
            ynew = cell(M, 1);
            s = 1;
            t = type;
        otherwise
            error('Invalid approx method');
    end

    st = 1;
    for j=s:s+M-1
        Xnew{j} = model.X{j};
        ynew{j} = model.y(st:st+size(Xnew{j}, 1)-1);
        st = st+size(Xnew{j}, 1);
    end
    
    Xnew{t} = [Xnew{t}; x];
    ynew{t} = [ynew{t}; y];
    
    