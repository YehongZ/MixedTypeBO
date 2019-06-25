function [ ret ] = getFuncValue(task, x, M, type)
    
    task_type = task{M+1};
    switch task_type
        case 'syn'
            seq = task{M+2};
            if numel(task) > M+2
                me = task{M+3};
                ret = task{type}(x, seq, me);
            else
                ret = task{type}(x, seq);
            end
        case 'Branin'
            ret = task{type}(x, 0);
        case 'Hartmann'
            me = task{M+3};
            ret = task{type}(x, me, 0);
        case 'LR'
%             x = transX(x, 'lr', true);
            ret = getFuncValue_lr(x, 'lr_model', type);
        case 'CNN'
%             x = transX(x, 'lr', true);
            ret = getFuncValue_cnn_bos(x, 'cnn_model', type);
        case 'Hawkes'
%             x = transX(x, 'lr', true);
            ret = getFuncValue_hawkes(x, 'haw_multi1', type);
        case 'RL_car'
            ret = getObsValue_car(x, type, task{M+2}, task{M+3});
        otherwise
            error('Invalid function type.');
    end