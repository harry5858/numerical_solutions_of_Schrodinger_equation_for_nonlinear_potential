clear;

%solve d2f/dt2+w(t)^2*f=0
%solve a system of equations
%rename f to y and t to x( not the space coordinate )
%system of equations
% dydx = z -------> g(x,y,z)
% dzdx = -1*w^2(x)*y----> h(x,y,z)

tmin=0;
tmax=12;
nt=120000;
dt=2*((tmax-tmin)/nt);
% t=zeros(1,nt);
% for i=1:nt
%     t(i)=i*dt;
% end
t=linspace(tmin,tmax,nt);
%  w initial and w final

w_i=2;
% w_f=6;
% 
% a=(w_f+w_i)/2;
% b=(w_f-w_i)/2;
% 
% dT=1;
w=2*ones(1,nt);
% w=a+(b*tanh((t-(tmax/2))/dT));
% figure(1)
% plot(t,w)

% RK 4
f=zeros(1,nt);
g=zeros(1,nt);

% initial condition
f(1)=1/sqrt(2*w_i);
g(1)=-1i*w_i*1/sqrt(2*w_i);
for n=1:2:nt-2
    k1=dt*(g(n));
    l1=dt*(-1*f(n)*w(n)^2);
    k2=dt*(g(n)+l1/2);
    l2=dt*(-1*(f(n)+k1/2)*w(n+1)^2);
    k3=dt*(g(n)+l2/2);
    l3=dt*(-1*(f(n)+k2/2)*w(n+1)^2);
    k4=dt*(g(n)+l3/2);
    l4=dt*(-1*(f(n)+k3/2)*w(n+2)^2);
    f(n+2)=f(n)+(1/6)*(k1+2*k2+2*k3+k4);
    g(n+2)=g(n)+(1/6)*(l1+2*l2+2*l3+l4);
end
nt1=nt/2;
t1=linspace(tmin,tmax,nt1);
f_new=zeros(1,nt1);
f_new(1)=f(1);
for i=2:nt1-1
f_new(i)=f(2*i+1);
end
f_new(nt1)=f_new(nt1-1);
figure(1)
plot(t1,real(f_new),'-b','Linewidth',1.5)
hold on
plot(t1,imag(f_new),'-r','Linewidth',1.5)
title('Real and Imaginary part of f(t) vs t')
xlabel('t')
set(gca,'fontsize',20)
legend('real part of f(t)','imag part of f(t)')
f_ex=(1/2)*exp(-1i*2*t1);
e1=real(f_new-f_ex);
figure(2)
plot(t1,e1)
% alpha=zeros(1,nt/2);
% beta=zeros(1,nt/2);
% for i=2:nt/2-1
% 
% f_out_a=sqrt(1/(2*w_f))*exp(-1i*w_f*t1(i));
% f_out_b=conj(f_out_a);
% df_out_a=-1i*w_f*f_out_a;
% df_out_b=1i*w_f*f_out_b;
% 
% df_in= (conj(f_new(i+1))-conj(f_new(i-1)))/(2*dt);
% 
% alpha(i)=1i*(conj(f_new(i))*df_out_a-df_in*f_out_a);
% beta(i)=1i*((conj(f_new(i))*df_out_b-df_in*f_out_b));
% 
% alpha(nt/2)=alpha(nt/2-1);
% alpha(1)=alpha(2);
% beta(nt/2)=beta(nt/2-1);
% beta(1)=beta(2);
% end
% 
% figure(2)
% 
% alpha_1=abs(alpha);
% beta_1=abs(beta);
% 
% plot(t1,alpha_1,'-b','Linewidth',1.5)
% hold on
% plot(t1,beta_1,'-r','Linewidth',1.5)
% title('absolute value of alpha and beta vs t')
% xlabel('t')
% set(gca,'fontsize',20)
% legend('abs(Alpha)','abs(Beta)')
% 
% figure(3)
% 
% alpha_real=real(alpha);
% alpha_imag=imag(alpha);
% 
% plot(t1,alpha_real,'Linewidth',1.5)
% hold on
% plot(t1,alpha_imag,'Linewidth',1.5)
% title('real and imaginary part of alpha vs t')
% xlabel('t')
% set(gca,'fontsize',20)
% 
% figure(5)
% 
% beta_real=real(beta);
% beta_imag=imag(beta);
% 
% plot(t1,beta_real,'Linewidth',1.5)
% hold on
% plot(t1,beta_imag,'Linewidth',1.5)
% title('real and imaginary part of beta vs t')
% xlabel('t')
% set(gca,'fontsize',20)
% 
% figure(6)
% 
% beta_2=beta_1.^2;
% 
% plot(t1,beta_2,'Linewidth',1.5)
% title('out Number operator (abs(beta))^2 vs t')
% xlabel('t')
% ylabel('out Number operator')
% set(gca,'fontsize',20)
% 
% g=zeros(1,nt/2);
% 
% for j=1:nt/2
% 
%    g(j)=abs(alpha(j))^2-abs(beta(j))^2;
%    
% end
% 
% figure(7)
% plot(t1,g,'Linewidth',1.5)
% title('abs(alpha)^2 - abs(beta)^2 vs t')
% xlabel('t')
% ylabel('abs(alpha)^2 - abs(beta)^2')
% set(gca,'fontsize',20)
% % Proability
% % Only consider the first 3 even terms for comparison with the Schrodinger
% % numerical solution
% % 
% alpha_p=alpha_1(nt/2-2000);
% beta_p=beta_1(nt/2-2000);
% 
% A=conj(beta_p)/(2*alpha_p);
% N=sqrt(1/abs(alpha_p));
% 
% C0=(abs(N))^2;
% C2=(abs(N*sqrt(2)*(A)))^2;
% C4=(abs(N*sqrt(6)*(A^2)))^2;
% 
% n=0:1:4;
% 
% figure(8)
% plot(n(1),C0,'b*','Markersize', 15)
% hold on
% plot(n(3),C2,'b*','Markersize', 15)
% plot(n(5),C4,'b*','Markersize', 15)
% axis([-0.1 5 0 1])
