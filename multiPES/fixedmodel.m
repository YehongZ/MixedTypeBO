function [ret] = fixedmodel(model, fixed)

    params = model.params;
    params(fixed.index) = fixed.value;
    ret = modelExpandParam(model, params);
    ret.params = params;