function [ ret ] = getObsValue(task, x, M, type)
    
    task_type = task{M+1};
    switch task_type
        case 'syn'
            seq = task{M+2};
            if numel(task) > M+2
                me = task{M+3};
                ret = task{type}(x, seq, me, true);
            else
                ret = task{type}(x, seq, true);
            end         
        case 'Branin'
            noise = task{M+2};
            ret = task{type}(x, noise(type));
        case 'Hartmann'
            noise = task{M+2};
            me = task{M+3};
            ret = task{type}(x, me, noise(type));
        otherwise
            error('Invalid function type.');
    end