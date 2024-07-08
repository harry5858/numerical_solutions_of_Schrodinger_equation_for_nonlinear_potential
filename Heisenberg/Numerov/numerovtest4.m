clear all;
%Numerov method
for nt=[100000 50000]
tmin=0;
tmax=10;
% nt=100000;
dt=(tmax-tmin)/nt;
t=linspace(tmin,tmax,nt);

w_i=2;
% w_f=6;
% a=(w_f+w_i)/2;
% b=(w_f-w_i)/2;
% dT=0.5;
% w=a+(b*tanh((t-(tmax/2))/dT));
w=2*ones(1,nt);

y=zeros(1,nt);

%initial condition
y(1)=1/sqrt(2*w_i);
dy(1)=-1i*w_i*1/sqrt(2*w_i);

%second initial value
a11=1+(dt^2)*(w(2)^2)/4;
a12=-1*(dt^2)*(w(3)^2)/24;
a21=-2+(5/6)*(dt^2)*(w(2)^2);
a22=1+(dt^2)*(w(3)^3)/12;
b1=y(1)+dt*dy(1)-(dt^2)*(7*(w(1)^2)*y(1))/24;
b2=-y(1)-(dt^2)*((w(1)^2)*y(1))/12;
y(2)=(a22*b1-a12*b2)/(a11*a22-a12*a21);
disp(y(2))

for n=2:nt-1

    y(n+1) = (2*(1-(5/12)*(dt^2)*(w(n)^2))*y(n)-(1+(1/12)*(dt^2)*(w(n-1)^2))*y(n-1))/(1+(1/12)*(dt^2)*(w(n+1)^2));
    
end
f_ex=(1/2)*exp(-1i*2*t);
% figure(1)
% plot(t,real(y),'r')
% hold on
% plot(t,imag(y),'b')
% plot(t,real(f_ex))

if nt==100000
e=zeros(1,nt);
for i=1:nt-1
e(i)=y(i+1)-f_ex(i);
end
% figure(2)
% plot(t,real(e))
% hold on
elseif nt==50000
e2=zeros(1,nt);
for i=1:nt-1
e2(i)=y(i+1)-f_ex(i);
end
end
end
t1=linspace(tmin,tmax,50000);
p=zeros(1,50000);
for i=1:50000
p=log2(abs(e2(i)/e(i+1)));
end
figure(100)
plot(t1,p)


































