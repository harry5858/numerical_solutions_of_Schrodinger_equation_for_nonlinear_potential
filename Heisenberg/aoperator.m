clear;

tmin=0;
tmax=12;
nt=120000;
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end

%  w initial and w final

w_i=3;
w_f=6;

a=(w_f+w_i)/2;
b=(w_f-w_i)/2;

dT=0.01;

w=a+b.*tanh((t-(tmax/2))/dT);

dwdt=1.5*(sech(t-6)).^2;

plot(t,w)

%solution array
u=zeros(1,nt);
v=zeros(1,nt);
%initial condition
u(1)=(w(2)+w(1))/(2*sqrt(w(2)*w(1)));
v(1)=(w(2)-w(1))/(2*sqrt(w(2)*w(1)));

for i=1:nt-2
    k1=dt*(-1i*w(i)*u(i)+(dwdt(i)/(2*w(i)))*conj(v(i)));
    l1=dt*(-1i*w(i)*v(i)+(dwdt(i)/(2*w(i)))*conj(u(i)));
    k2=dt*(-1i*w(i+1)*(u(i)+k1/2)+(dwdt(i+1)/(2*w(i+1)))*conj(v(i)+l1/2));
    l2=dt*(-1i*w(i+1)*(v(i)+k1/2)+(dwdt(i+1)/(2*w(i+1)))*conj(u(i)+l1/2));
    k3=dt*(-1i*w(i+1)*(u(i)+k2/2)+(dwdt(i+1)/(2*w(i+1)))*conj(v(i)+l2/2));
    l3=dt*(-1i*w(i+1)*(v(i)+k2/2)+(dwdt(i+1)/(2*w(i+1)))*conj(u(i)+l2/2));
    k4=dt*(-1i*w(i+2)*(u(i)+k3/2)+(dwdt(i+2)/(2*w(i+2)))*conj(v(i)+l3/2));
    l4=dt*(-1i*w(i+2)*(v(i)+k3/2)+(dwdt(i+2)/(2*w(i+2)))*conj(u(i)+l3/2));
    u(i+1)=u(i)+(1/6)*(k1+2*k2+2*k3+k4);
    v(i+1)=v(i)+(1/6)*(l1+2*l2+2*l3+l4);
    u(nt)=u(nt-1);
    v(nt)=v(nt-1);
end
u_real=real(u);
u_imag=imag(u);
v_real=real(v);
v_imag=imag(v);
figure(1)
plot(t,u_real)
hold on
plot(t,u_imag)
plot(t,abs(u))

figure(2)
plot(t,v_real)
hold on
plot(t,v_imag)
plot(t,abs(v))

norm_check=zeros(1,nt);
for j=1:nt

   norm_check(j)=abs(u(j))^2-abs(v(j))^2;
   
end
figure(3)
plot(t,norm_check)

% Proability
% Only consider the first 3 even terms for comparison with the Schrodinger
% numerical solution

alpha_p=u(nt);
beta_p=v(nt);

A=-1*conj(beta_p)/(2*alpha_p);
N=sqrt(1/abs(alpha_p));

C0=(abs(N))^2;
C2=(abs(N*sqrt(2)*(A)))^2;
C4=(abs(N*sqrt(6)*(A^2)))^2;

n=0:1:4;

figure(4)
plot(n(1),C0,'b*','Markersize', 15)
hold on
plot(n(3),C2,'b*','Markersize', 15)
plot(n(5),C4,'b*','Markersize', 15)
axis([-0.1 5 0 1])
xlabel('n','Fontsize', 24)
ylabel('probability |Cn|^2','Fontsize', 24)
title('probability |<2n out | 0 in>|^2 (sudden transition) ','Fontsize', 20)

