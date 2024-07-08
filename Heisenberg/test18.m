clear all;
%solve d2f/dt2+w(t)^2*f=0
%solve a system of equations
%rename f to y and t to x( not the space coordinate )
%system of equations
% dydx = z -------> g(x,y,z)
% dzdx = w^2(x)*y----> h(x,y,z)

tmin=0;
tmax=12;
nt=600000;
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end

%  w initial and w final

w_i=3;
w_f=1:0.5:20;
n_i=numel(w_f);
wfwi=zeros(1,n_i);

for n=1:n_i
wfwi(n)=w_f(n)/w_i;
end

% C0=zeros(1,n_i);
C2=zeros(1,n_i);
% C4=zeros(1,n_i);

for j=1:n_i
    
a=(w_f(j)+w_i)/2;
b=(w_f(j)-w_i)/2;
dT=1/2;
% RK 4

f=zeros(1,nt);
g=zeros(1,nt);

% initial condition
f(1) = sqrt(1/6);
g(1) = -1i*sqrt(3/2);

for n=1:nt-2
    
    t(n+1) = t(n)+dt;
    k1=dt*(g(n));
    l1=dt*(-1*f(n)*(a+b*tanh((t(n)-6)/dT))^2);
    k2=dt*(g(n)+l1/2);
    l2=dt*(-1*(f(n)+k1/2)*(a+b*tanh((t(n+1)-6)/dT))^2);
    k3=dt*(g(n)+l2/2);
    l3=dt*(-1*(f(n)+k2/2)*(a+b*tanh((t(n+1)-6)/dT))^2);
    k4=dt*(g(n)+l3/2); 
    l4=dt*(-1*(f(n)+k3/2)*(a+b*tanh((t(n+2)-6)/dT))^2);
    f(n+1)=f(n)+(1/6)*(k1+2*k2+2*k3+k4);
    g(n+1)=g(n)+(1/6)*(l1+2*l2+2*l3+l4);
    f(nt)=f(nt-1);
    g(nt)=g(nt-1);
end

% figure(1)
% plot(x,real(y),'-b','Linewidth',1.5)
% hold on
% plot(x,imag(y),'-r','Linewidth',1.5)
% xlabel('t')
% title('real and imaginary part of f(t) vs t')
% set(gca,'fontsize',20)
% % legend('real part of f(t)','imag part of f(t)')
% hold off
alpha=zeros(1,nt);
beta=zeros(1,nt);

for i=2:nt-1

f_out_a=sqrt(1/(2*w_f(j)))*exp(-1i*w_f(j)*t(i));
f_out_b=conj(f_out_a);
df_out_a=-1i*w_f(j)*f_out_a;
df_out_b=1i*w_f(j)*f_out_b;

df_in= (conj(f(i+1))-conj(f(i-1)))/(2*dt);

alpha(i)=1i*(conj(f(i))*df_out_a-df_in*f_out_a);
beta(i)=1i*((conj(f(i))*df_out_b-df_in*f_out_b));

alpha(nt-1)=alpha(nt-2);
alpha(nt)=alpha(nt-1);
alpha(1)=alpha(2);
beta(nt-1)=beta(nt-2);
beta(nt)=beta(nt-1);
beta(1)=beta(2);
end

% figure(2)

% alpha_1=abs(alpha);

% plot(x,alpha_1,'Linewidth',1.5)
% title('absolute value of alpha vs t')
% xlabel('t')
% ylabel('abs(alpha)')
% set(gca,'fontsize',20)

% figure(3)

% alpha_real=real(alpha);
% alpha_imag=imag(alpha);

% plot(x,alpha_real,'Linewidth',1.5)
% hold on
% plot(x,alpha_imag,'Linewidth',1.5)
% title('real and imaginary part of alpha vs t')
% xlabel('t')
% set(gca,'fontsize',20)

% figure(4)

% beta_1=abs(beta);

% plot(x,beta_1,'Linewidth',1.5)
% title('absolute vlaue of beta vs t')
% xlabel('t')
% ylabel('abs(beta)')
% set(gca,'fontsize',20)

% figure(5)

% beta_real=real(beta);
% beta_imag=imag(beta);

% plot(x,beta_real,'Linewidth',1.5)
% hold on
% plot(x,beta_imag,'Linewidth',1.5)
% title('real and imaginary part of beta vs t')
% xlabel('t')
% set(gca,'fontsize',20)

% figure(6)

% beta_2=beta_1.^2;

% plot(x,beta_2,'Linewidth',1.5)
% title('out Number operator (abs(beta))^2 vs t')
% xlabel('t')
% ylabel('out Number operator')
% set(gca,'fontsize',20)

% g=zeros(1,nx);
% 
% for m=1:nx
% 
%    g(m)=abs(alpha(m))^2-abs(beta(m))^2;
%    
% end

% figure(7)

% plot(x,g,'Linewidth',1.5)
% hold on
% title('abs(alpha)^2 - abs(beta)^2 vs t')
% xlabel('t')
% ylabel('abs(alpha)^2 - abs(beta)^2')
% set(gca,'fontsize',20)

% Proability
% Only consider the first 3 even terms for comparison with the Schrodinger
% numerical solution

alpha_p=alpha(nt-10000);
beta_p=beta(nt-10000);

A=-1*conj(beta_p)/(2*alpha_p);
N=sqrt(1/abs(alpha_p));

% C0(j)=(abs(N))^2;
C2(j)=(abs(N*sqrt(2)*(A)))^2;
% C4(j)=(abs(N*sqrt(6)*(A^2)))^2;

end
figure(8)
% semilogy(dw,C0,'b*','Markersize', 15)
% hold on
semilogy(wfwi,C2,'r')
% semilogy(dw,C4,'k')
% axis([-3 7.5 0 1])
% legend('C0','C2','C4')
xlabel('w_f/w_i','Fontsize', 24)
ylabel('probability |Cn|^2','Fontsize', 24)
% title('probability |<2n out | 0 in>|^2 , n=0,1,2 ','Fontsize', 20)
