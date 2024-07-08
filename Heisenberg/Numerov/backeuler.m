clear all;
for nt=[200000 100000]
tmin=0;
tmax=100;
% nt=10000;
t=linspace(tmin,tmax,nt);
dt=tmax/nt;
% set up for solution and predictor solution
y1=zeros(1,nt);
y12=zeros(1,nt);
y2=zeros(1,nt);
y22=zeros(1,nt);
w_i=2;
% w_f=2;
% a=(w_f+w_i)/2;
% b=(w_f-w_i)/2;
% dT=0.5;
% w=a+(b*tanh((t-(tmax/2))/dT));
w=w_i*ones(1,nt);
%initial condition
y1(1)=1/sqrt(2*w_i);
y2(1)=-1i*w_i*1/sqrt(2*w_i);
y12(1)=1/sqrt(2*w_i);
y22(1)=-1i*w_i*1/sqrt(2*w_i);
% Predictor Forward Euler
% y12(2)=y1(1)+dt*y2(1);
% y22(2)=y2(2)-dt*(w(1)^2)*(y1(1));
for i=1:nt-1
    y12(i+1)=y12(i)+dt*y22(i)-0.5*(dt^2)*(w(i)^2)*y12(i);
    y22(i+1)=y22(i)-dt*(w(i)^2)*(y12(i));
end
% Corrector Backward Euler
% y1(2)=y1(1)+dt*y22(2);
% y2(2)=y2(1)-(w(2)^2)*y12(2);

% for n=1:nt-1
%     y1(n+1)=y1(n)+dt*y22(n+1);
%     y2(n+1)=y2(n)-dt*(w(n+1)^2)*y12(n+1);
% end
% 
% figure(1)
% plot(t,real(y1))
% hold on
% plot(t,imag(y1))
% f_ex=0.5*exp(-1i*2*t);
% e=y1-f_ex;
% figure(2)
% plot(t,e)

% Corrector AM2
% y1(2)=y1(1)+0.5*dt*(y22(2)+y2(1));
% y2(2)=y2(1)-0.5*dt*((w(2)^2)*y12(2)+(w(1)^2)*y1(1));

for n=1:nt-1
    y1(n+1)=y1(n)+0.5*dt*(y22(n+1)+y2(n));
    y2(n+1)=y2(n)-0.5*dt*((w(n+1)^2)*y12(n+1)+(w(n)^2)*y1(n));
end

% figure(1)
% plot(t,real(y1))
% hold on
% plot(t,imag(y1))

f_ex=0.5*exp(-1i*2*t);
%     if nt==200000
    e=zeros(1,nt);
    for i=1:nt-1
    e(i)=y1(i+1)-f_ex(i);
    end
    figure(2)
    plot(t,real(e))
    hold on
    drawnow;
%     elseif nt==100000
%         e2=zeros(1,nt);
%     for i=1:nt-1
%     e2(i)=y1(i+1)-f_ex(i);
%     end
%     end
end
% t1=linspace(tmin,tmax,100000);
% p=zeros(1,100000);
% for i=1:100000
%     p(i)=log2(abs(e2(i)/e(i+1)));
% end
% figure(3)
% plot(t1,p)









