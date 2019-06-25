function [ret] = get_cost(x, type)
	epoch = x(1);
	if type == 1
		ret = 1;
	else
		ret = 1/(5*epoch/20);
	end
end
