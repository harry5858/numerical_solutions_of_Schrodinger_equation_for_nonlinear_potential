function xdot=myfunc(t,x,w_f,w_i,dT)
xdot = [x(2); -((((w_f+w_i)/2)+((w_f-w_i)/2).*tanh((t-6)/dT)).^2).*x(1)];
end
