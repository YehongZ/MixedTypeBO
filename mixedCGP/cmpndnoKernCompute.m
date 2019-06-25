function k = cmpndnoKernCompute(kern, x, x2)
% Compute the cmpnd kernel without adding the noise variance.
% Input:
%   kern: the kernel structure for which the matrix is computed.
%   x: the input matrix associated with the rows of the kernel.
%   x2: the input matrix associated with the columns of the kernel.
% Output:
%	  k - the kernel matrix computed at the given inputs.
%
% Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
% This function is modified from CMPNDKERNCOMPUTE by Yehong	

if nargin > 2
  i = 1;
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved in the kernel.
    k = kernCompute(kern.comp{i}, ...
                         x(:, kern.comp{i}.index), ...
                         x2(:, kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    k = kernCompute(kern.comp{i}, x, x2);
  end
  for i = 2:length(kern.comp)-1     %%%%% exclude the white noise
    if ~isempty(kern.comp{i}.index)
      % only part of the data is involved in the kernel.
      k  = k + kernCompute(kern.comp{i}, ...
                           x(:, kern.comp{i}.index), ...
                           x2(:, kern.comp{i}.index));
    else
      % all the data is involved with the kernel.
      k  = k + kernCompute(kern.comp{i}, x, x2);
    end
  end
else
  i = 1;
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved with the kernel.
    k  = kernCompute(kern.comp{i}, x(:, kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    k  = kernCompute(kern.comp{i}, x);
  end
  for i = 2:length(kern.comp)-1     %%%%% exclude the white noise
    if ~isempty(kern.comp{i}.index)
      % only part of the data is involved with the kernel.
      k  = k + kernCompute(kern.comp{i}, x(:, kern.comp{i}.index));
    else
      % all the data is involved with the kernel.
      k  = k + kernCompute(kern.comp{i}, x);
    end
  end
end
