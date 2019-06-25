function [Xnew, ynew, noise] = updateMixedXY(model, x, y, type)

    nlf = model.nlf;
    
    switch model.approx
        case 'ftc'
            Xnew = model.X;
            ynew = model.oY;
            noise = model.noise;
            t = type + nlf;
        case {'dtc','fitc','pitc'}
            error('Not implement yet');
        otherwise
            error('Invalid approx method');
    end
    
    Xnew{t} = [Xnew{t}; x];
    ynew{t} = [ynew{t}; y];
    noise{type} = [noise{type}; 0];
    