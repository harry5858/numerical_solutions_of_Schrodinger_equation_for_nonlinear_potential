clear all;
%leapfrog
for nt=[100000 50000]
tmin=0;
tmax=10;
% nt=10000;
t=linspace(0,10,nt);
dt=(tmax-tmin)/nt;

w_i=2;
% w_f=6;
% dT=0.5;
% a=(w_f+w_i)/2;
% b=(w_f-w_i)/2;
% w=a+b.*tanh((t-tmax/2)/dT);
w=2*ones(1,nt);
y1=zeros(1,nt);
y2=zeros(1,nt);

%initial condition

y1(1)=sqrt(1/(2*w_i));
y2(1)=-1i*sqrt(w_i/2);

for n=1:nt-1

y1(n+1)=y1(n)+dt*y2(n)-0.5*(dt^2)*((w(n)^2)*y1(n));

y2(n+1)=y2(n)-0.5*dt*((w(n)^2)*y1(n)+(w(n+1)^2)*y1(n+1));

end

figure(123)
plot(t,real(y1),'b')
hold on
plot(t,imag(y1),'r')
f_ex=0.5*exp(-1i*2*t);

% plot(t,real(f_ex),'k')

    if nt==100000
    e=zeros(1,nt);
    for i=1:nt-1
    e(i)=y1(i+1)-f_ex(i);
    end
%     figure(2)
%     plot(t,real(e))
%     hold on
%     drawnow;
    elseif nt==50000
        e2=zeros(1,nt);
    for i=1:nt-1
    e2(i)=y1(i+1)-f_ex(i);
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







