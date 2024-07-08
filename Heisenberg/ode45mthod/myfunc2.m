function dydt=myfunc2(t,y,w_i,w_f,dT)
a=(w_f+w_i)/2;
b=(w_f-w_i)/2;
w=a+(b*tanh((t-(6))/dT));

dydt = [y(2);-(w^2)*y(1)];
