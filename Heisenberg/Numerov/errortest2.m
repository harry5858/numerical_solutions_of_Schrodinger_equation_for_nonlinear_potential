%Numerov method
for nt=[100000 50000 25000]
    
    tmin=0;
    tmax=10;
    dt=(tmax-tmin)/nt;
    t=zeros(1,nt);
    for i=1:nt
    t(i)=i*dt;
    end

    w_i=2;
    w_f=2;

    a=(w_f+w_i)/2;
    b=(w_f-w_i)/2;

    dT=1;

    w=2*ones(1,nt);

    f=zeros(1,nt);

%initial condition
f(1)=1/sqrt(2*w_i);
g(1)=-1i*w_i*1/sqrt(2*w_i);

%RK 2 for another initial condition

k1=dt*(g(1));
l1=dt*(-1*f(1)*w(1)^2);
k2=dt*(g(1)+l1/2);
l2=dt*(-1*(f(1)+k1/2)*w(2)^2);
k3=dt*(g(1)+l2/2);
l3=dt*(-1*(f(1)+k2/2)*w(2)^2);
k4=dt*(g(1)+l3/2);
l4=dt*(-1*(f(1)+k3/2)*w(3)^2);
f1(3)=f(1)+(1/6)*(k1+2*k2+2*k3+k4);
f(2)=f1(3);
% 
%     t(2) = t(1)+dt;
%     f(2)=f(1)+dt*g(1)-0.5*(dt^2)*(w(1)^2)*f(1);

% new set of initial conditions are f(1) and f(2)

for n=2:nt-1
    
f(n+1) = (2*(1-(5/12)*(dt^2)*(w(n)^2))*f(n)-(1+(1/12)*(dt^2)*(w(n-1)^2))*f(n-1))/(1+(1/12)*(dt^2)*(w(n+1)^2));

end

f(nt)=f(nt-1);

% plot(t,real(f))
% hold on
% plot(t,imag(f))
f_ex=(1/2)*exp(-1i*2*t);

if nt==100000
    f1=f;
    f1_ex=f_ex;
elseif nt==50000
    f2=f;
    f2_ex=f_ex;
elseif nt==25000
    f4=f;
    f4_ex=f_ex;
end
end
t1=linspace(0,12,100000);
t2=linspace(0,12,50000);
t4=linspace(0,12,25000);

e1=zeros(1,100000);
e2=zeros(1,50000);
e4=zeros(1,25000);
for i=1:100000-1
e1(i)=f1(i+1)-f1_ex(i);
end
for i=1:50000-1
e2(i)=f2(i+1)-f2_ex(i);
end
for i=1:25000-1
e4(i)=f4(i+1)-f4_ex(i);
end
figure(1)
plot(t1,e1)
hold on
plot(t2,e2)

plot(t4,e4)

xlabel('t')
ylabel('error')
title('error vs t')
set(gca,'fontsize',20)



