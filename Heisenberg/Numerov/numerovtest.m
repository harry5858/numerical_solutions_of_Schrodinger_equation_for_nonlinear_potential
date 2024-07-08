%Numerov method
clear all;
tmin=0;
tmax=12;
nt=240000;
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end

w_i=3;
w_f=6;

a=(w_f+w_i)/2;
b=(w_f-w_i)/2;

dT=0.0001;

w=a+(b*tanh((t-(tmax/2))/dT));

%initial condition
f(1)=1/sqrt(2*w_i);
g(1)=-1i*w_i*1/sqrt(2*w_i);

%get new initial value
f(2)=f(1)+dt*g(1)-0.5*(dt^2)*(w(1)^2)*f(1);

% new set of initial conditions is f(1) and f(2)

for n=2:nt-1
    
f(n+1) = (2*(1-(5/12)*(dt^2)*(w(n)^2))*f(n)-(1+(1/12)*(dt^2)*(w(n-1)^2))*f(n-1))/(1+(1/12)*(dt^2)*(w(n+1)^2));

end

f(nt)=f(nt-1);
figure(1)
plot(t,real(f),'b')
hold on
plot(t,imag(f),'r')

alpha=zeros(1,nt);
beta=zeros(1,nt);
for i=2:nt-1

f_out_a=sqrt(1/(2*w_f))*exp(-1i*w_f*t(i));
f_out_b=conj(f_out_a);
df_out_a=-1i*w_f*f_out_a;
df_out_b=1i*w_f*f_out_b;

df_in= ((f(i+1))-(f(i-1)))/(2*dt);

alpha(i)=-1i*(f(i)*df_out_b-df_in*f_out_b);
beta(i)=1i*(f(i)*df_out_a-df_in*f_out_a);

alpha(nt-1)=alpha(nt-2);
alpha(nt)=alpha(nt-1);
alpha(1)=alpha(2);
beta(nt-1)=beta(nt-2);
beta(nt)=beta(nt-1);
beta(1)=beta(2);
end

figure(2)

alpha_1=abs(alpha);
beta_1=abs(beta);

plot(t,alpha_1,'-b','Linewidth',1.5)
hold on
plot(t,beta_1,'-r','Linewidth',1.5)
title('absolute value of alpha and beta vs t')
xlabel('t')
set(gca,'fontsize',20)
legend('abs(Alpha)','abs(Beta)')

figure(3)

alpha_real=real(alpha);
alpha_imag=imag(alpha);

plot(t,alpha_real,'Linewidth',1.5)
hold on
plot(t,alpha_imag,'Linewidth',1.5)
title('real and imaginary part of alpha vs t')
xlabel('t')
set(gca,'fontsize',20)

figure(5)

beta_real=real(beta);
beta_imag=imag(beta);

plot(t,beta_real,'Linewidth',1.5)
hold on
plot(t,beta_imag,'Linewidth',1.5)
title('real and imaginary part of beta vs t')
xlabel('t')
set(gca,'fontsize',20)

figure(6)

beta_2=beta_1.^2;

plot(t,beta_2,'Linewidth',1.5)
title('out Number operator (abs(beta))^2 vs t')
xlabel('t')
ylabel('out Number operator')
set(gca,'fontsize',20)

g=zeros(1,nt);

for j=1:nt

   g(j)=abs(alpha(j))^2-abs(beta(j))^2;
   
end

% figure(7)
plot(t,g,'Linewidth',1.5)
title('abs(alpha)^2 - abs(beta)^2 vs t')
xlabel('t')
ylabel('abs(alpha)^2 - abs(beta)^2')
set(gca,'fontsize',20)
% Proability
% Only consider the first 3 even terms for comparison with the Schrodinger
% numerical solution

alpha_p=alpha_1(nt-3);
beta_p=beta_1(nt-3);

A=conj(beta_p)/(2*alpha_p);
N=sqrt(1/abs(alpha_p));

C0=(abs(N))^2;
C2=(abs(N*sqrt(2)*(A)))^2;
C4=(abs(N*sqrt(6)*(A^2)))^2;

n=0:1:4;

figure(8)
plot(n(1),C0,'bo','Markersize', 15)
hold on
plot(n(3),C2,'bo','Markersize', 15)
plot(n(5),C4,'bo','Markersize', 15)
txt0 = '\leftarrow |C_0|^2 = 0.9428';
txt2 = '\leftarrow |C_2|^2 = 0.05237';
txt4 = '\downarrow |C_4|^2 = 0.004364';
t0=text(0,C0,txt0);
t2=text(2,C2,txt2);
t4=text(4,C4+0.03,txt4);
t0.FontSize = 20;
t2.FontSize = 20;
t4.FontSize = 20;
xlabel('n','Fontsize', 24)
title('Transition probabilities |C(n)|^2 for n=0,1,2,3,4 ')
axis([-0.1 5 0 1])
set(gca,'fontsize',20)


