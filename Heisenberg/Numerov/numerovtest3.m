%Numerov method
clear all;
tmin=0;
tmax=10;
nt=200000;
nt1=nt/2;
t=linspace(tmin,tmax,nt);
t1=linspace(tmin,tmax,nt1);
dt=2*tmax/nt;
dt2=tmax/nt;

w_i=3;
w_f=6;
a=(w_f+w_i)/2;
b=(w_f-w_i)/2;
dT=1;
w=a+(b*tanh((t-(tmax/2))/dT));

%initial condition
% f(1)=1/sqrt(2*w_i);
% g(1)=-1i*w_i*1/sqrt(2*w_i);

f=zeros(1,nt);

% initial condition
y_d(1)=1/sqrt(2*w_i);
y1_d(1)=-1i*w_i*1/sqrt(2*w_i);

k0= (dt^2)*(-1*(w(1)^2)*y_d(1));
k1= (dt^2)*(-1*(w(2)^2)*(y_d(1)+dt2*y1_d(1)+k0/8));
k2= (dt^2)*(-1*(w(3)^2)*(y_d(1)+dt*y1_d(1)+0.5*k1));

y_d(3)=y_d(1)+dt*y1_d(1)+(k0+2*k1)/6;

%get new initial value
% f(2)=f(1)+dt*g(1)-0.5*(dt^2)*(w(1)^2)*f(1);

% new set of initial conditions is f(1) and f(2)
f(1)=1/sqrt(2*w_i);
f(2)=y_d(3);

for n=2:nt-1
    
% f(n+1) = (2*(1-(5/12)*((dt2)^2)*(w(n)^2))*f(n)-(1+(1/12)*((dt2)^2)*(w(n-1)^2))*f(n-1))/(1+(1/12)*((dt2)^2)*(w(n+1)^2));

fn1=2*f(n)*(1-(5/12)*(dt2)^2*(w(n)^2));
fn2=f(n-1)*(1+(1/12)*(dt2)^2*(w(n-1)^2));
fd=(1+(1/12)*(dt2)^2*(w(n-1)^2));
f(n+1)= (fn1-fn2)/fd;

end

f(nt)=f(nt-1);

plot(t,real(f))
hold on
plot(t,imag(f))