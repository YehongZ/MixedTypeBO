function ynew = transY(y, ymean, method, inverse)
    
    if nargin < 4
        inverse = false;
    end
    
    if inverse
        switch method
            case 'lr'         
                ynew = y + ymean;
            case 'cnn'
                ynew = y + ymean;
            otherwise
                ynew = y;
        end
    else
        switch method
            case 'lr'         
                ynew = y - ymean;
            case 'cnn'
                ynew = y - ymean;
            otherwise
                ynew = y;
        end
    end