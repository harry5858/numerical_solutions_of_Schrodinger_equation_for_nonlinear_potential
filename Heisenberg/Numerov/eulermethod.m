clear all;
for nt=[10000 5000]
tmin=0;
tmax=10;
% nt=10000;
t=linspace(tmin,tmax,nt);
dt=tmax/nt;
y1=zeros(1,nt);
y2=zeros(1,nt);
w_i=2;
w_f=2;
a=(w_f+w_i)/2;
b=(w_f-w_i)/2;
dT=0.5;
w=a+(b*tanh((t-(tmax/2))/dT));
% w=2*ones(1,nt);
dw=(b/dT).*(sech((t-(tmax/2))./dT)).^2;
% ddw=(-2*b/(dT^2))*((sech*((t-(tmax/2))/dT))^2)*(tanh((t-(tmax/2))/dT));
%initial condition
y1(1)=1/sqrt(2*w_i);
y2(1)=-1i*w_i*1/sqrt(2*w_i);
for n=1:nt-1
y1(n+1)=y1(n)+dt*y2(n)-0.5*(dt^2)*(w(n)^2)*y1(n)-(1/6)*(dt^3)*(2*w(n)*dw(n)+(w(n)^2)*y2(n));
y2(n+1)=y2(n)-dt*(w(n)^2)*y1(n)-0.5*(dt^2)*(2*w(n)*dw(n)+(w(n)^2)*y2(n));               %+(1/6)*(dt^3)*(w(n)^4)*y1(n);
end

% figure(1)
% plot(t,real(y1))
% hold on
% plot(t,imag(y1))
f_ex=0.5*exp(-1i*2*t);
    if nt==10000
    e=zeros(1,nt);
    for i=1:nt-1
    e(i)=y1(i+1)-f_ex(i);
    end
%     figure(2)
%     plot(t,real(e))
%     hold on
%     drawnow;
    elseif nt==5000
        e2=zeros(1,nt);
    for i=1:nt-1
    e2(i)=y1(i+1)-f_ex(i);
    end
    end
end
t1=linspace(tmin,tmax,5000);
p=zeros(1,5000);
for i=1:5000
    p(i)=log2(abs(e2(i)/e(i+1)));
end
figure(3)
plot(t1,p)













