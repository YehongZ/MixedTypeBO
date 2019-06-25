% we use the predictive mean of a trained model to simulate the true function
% value
function [ ret ] = getFuncValue_cnn_bos(x, name, t)

    load(name)
    x = transX(x, 'cnn');
    nlf = model.nlf;
    [mu, ~] = multigpMixedPosteriorMeanVar(model, x);
    ret = mu{t+nlf};
    