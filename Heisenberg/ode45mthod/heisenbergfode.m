%use ode45 to solve the equation of motion
clear all;
trange=linspace(0,12,120000);
% dt=max(trange)/60000;
% nt=numel(trange);
w_i=2;
w_f=2;
dT=0.0001;
[t, y] = ode45(@(t,y) myfunc2(t,y,w_i,w_f,dT),trange,[sqrt(1/(2*w_i));-1i*(w_i)*sqrt(1/(2*w_i))]);
figure(10)
plot(t,real(y(:,1)),'b')
hold on
plot(t,imag(y(:,1)),'r')

y_new=y(:,1).';
y_ex=(1/2)*exp(-1i*2*trange);
% e1=real(y_new-y_ex);
figure(2)
plot(trange,real(y_ex))