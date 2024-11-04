function [newx]=RemoveLinear(n,x,t)

%     Roberto Reno', September 2003

%     input: the function x_i(t_i) i=1,n
%     output newx_i such that x_1 = x_n (removal of a linear trend)
      
      %for i=1:n
      %   newx(i) = x(i) - (x(n)-x(1))/(t(n)-t(1))*(t(i)-t(1));
      %end
      
      newx = x - (x(n)-x(1))/(t(n)-t(1))*(t-t(1));