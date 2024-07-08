%Numerov method
clear all;
for nt=[100000 50000]
    tmin=0;
    tmax=100;
    dt=(tmax-tmin)/nt;
    t=zeros(1,nt);
    for i=1:nt
        t(i)=i*dt;
    end

    w_i=2;
    % w_f=6;
    % 
    % a=(w_f+w_i)/2;
    % b=(w_f-w_i)/2;
    % 
    % dT=0.5;
    % 
    % w=a+(b*tanh((t-(tmax/2))/dT));
    w=w_i*ones(1,nt);
    f=zeros(1,nt);

    %initial condition
    f(1)=1/sqrt(2*w_i);
    g(1)=-1i*w_i*1/sqrt(2*w_i);

    %RK 2 for another initial condition

    t(2) = t(1)+dt;
    f(2)=f(1)+dt*g(1)-0.5*(dt^2)*(w(1)^2)*f(1);
        
%     a11=1+(dt^2)*(w(2)^2)/4;
%     a12=-1*(dt^2)*(w(3)^2)/24;
%     a21=-2+(5/6)*(dt^2)*(w(2)^2);
%     a22=1+(dt^2)*(w(3)^3)/12;
%     b1=f(1)+dt*g(1)-(dt^2)*(7*(w(1)^2)*f(1))/24;
%     b2=-f(1)-(dt^2)*((w(1)^2)*f(1))/12;
%     f(2)=(a22*b1-a12*b2)/(a11*a22-a12*a21);

    % new set of initial conditions are f(1) and f(2)

    for n=2:nt-1

    f(n+1) = (2*(1-(5/12)*(dt^2)*(w(n)^2))*f(n)-(1+(1/12)*(dt^2)*(w(n-1)^2))*f(n-1))/(1+(1/12)*(dt^2)*(w(n+1)^2));

    end

    f_ex=(1/2)*exp(-1i*2*t);

    % figure(1)
    % plot(t,real(f))
    % hold on
    % plot(t,imag(f))

    if nt==100000
    e=zeros(1,nt);
    for i=1:nt-1
    e(i)=f(i+1)-f_ex(i);
    end
%     figure(2)
%     plot(t,real(e))
%     xlabel('t')
%     ylabel('error')
%     title('Error vs t')
%     set(gca,'fontsize',20)
%     legend('dt=0.005','dt=0.01')
%     hold on
%     drawnow;
    elseif nt==50000
        e2=zeros(1,nt);
    for i=1:nt-1
    e2(i)=f(i+1)-f_ex(i);
    end
    end
end
t1=linspace(tmin,tmax,50000);
p=zeros(1,50000);
for i=1:50000
    p(i)=log2(abs(e2(i)/e(i+1)));
end
figure(3)
plot(t1,p)
xlabel('t')
ylabel('p')
title('Order of accuracy vs t')
set(gca,'fontsize',20)





