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
w_f=0:0.5:10;
n_i=numel(w_f);
dw=zeros(1,n_i);

for n=1:n_i
dw(n)=w_f(n)-w_i;
end

% C0=zeros(1,n_i);
C2=zeros(1,n_i);
% C4=zeros(1,n_i);

for k=1:n_i
    
a=(w_f(k)+w_i)/2;
b=(w_f(k)-w_i)/2;

dT=1;

w=a+b.*tanh((t-(tmax/2))/dT);

% dwdt=1.5*(sech(t-6)).^2;
% 
% plot(t,w)

A=zeros(1,nt);
B=zeros(1,nt);
u=zeros(1,nt);
v=zeros(1,nt);
ut=zeros(1,nt);
vt=zeros(1,nt);

for i=1:nt-1
    A(i)=(w(i+1)+w(i))/(2*sqrt(w(i)*w(i+1)));
    B(i)=(w(i+1)-w(i))/(2*sqrt(w(i)*w(i+1)));
    A(nt)=A(nt-1);
    B(nt)=B(nt-1);
end

u(1)=A(1);
v(1)=B(1);

for j=2:nt-1
   u(j)=(A(j)*u(j-1)+B(j)*conj(v(j-1)))*cos(w(j)*(t(j)-t(j-1)))-1i*(A(j)*u(j-1)-B(j)*conj(v(j-1)))*sin(w(j)*(t(j)-t(j-1)));
   v(j)=(A(j)*v(j-1)+B(j)*conj(u(j-1)))*cos(w(j)*(t(j)-t(j-1)))-1i*(A(j)*v(j-1)-B(j)*conj(u(j-1)))*sin(w(j)*(t(j)-t(j-1)));
   u(nt)=u(nt-1);
   v(nt)=v(nt-1);
end

for m=2:nt-1
    ut(m)=u(m)*exp(-1i*w(m+1)*(t(m)-t(m-1)));
    vt(m)=v(m)*exp(1i*w(m+1)*(t(m)-t(m-1)));
    ut(nt)=ut(nt-1);
    vt(nt)=vt(nt-1);
end


% figure(1)
% plot(t,real(ut))
% hold on
% plot(t,imag(ut))
% figure(2)
% plot(t,real(vt))
% hold on
% plot(t,imag(vt))
% figure(3)
% plot(t,abs(ut))
% figure(4)
% plot(t,abs(vt))
% figure(5)
% plot(t,(abs(ut)).^2-(abs(vt)).^2)

% Proability
% Only consider the first 3 even terms for comparison with the Schrodinger
% numerical solution

alpha_p=ut(nt);
beta_p=vt(nt);

A=-1*conj(beta_p)/(2*alpha_p);
N=sqrt(1/abs(alpha_p));

% C0=(abs(N))^2;
C2(k)=(abs(N*sqrt(2)*(A)))^2;
% C4=(abs(N*sqrt(6)*(A^2)))^2;
% 

end

% n=0:1:4;
figure(6)
% plot(n(1),C0,'b*','Markersize', 15)
% hold on
semilogy(dw,C2,'b')
% plot(n(5),C4,'b*','Markersize', 15)
xlabel('n','Fontsize', 24)
ylabel('probability |Cn|^2','Fontsize', 24)
title('probability |<2n out | 0 in>|^2 (sudden transition) ','Fontsize', 20)
