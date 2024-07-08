clear all;
%solve d2f/dt2+w(t)^2*f=0
%solve a system of equations
%rename f to y and t to x( not the space coordinate )
%system of equations
% dydx = z -------> g(x,y,z)
% dzdx = w^2(x)*y----> h(x,y,z)

tmin=0;
tmax=24;
nt=120000;
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end

%  w initial and w final

c_i=2;
c_f=4;
a=(c_f+c_i)/2;
b=(c_f-c_i)/2;
%c=a+b*tanh(t-6);
k=0:0.05:5;
% for i=1:11
% dT=1.5;
% w2=k(i)^2+(a+b.*tanh((t-6)/dT)).^2;
% hold on
% figure(1)
% plot(t,w2)
% drawnow;
% end
n_i=numel(k);
% alpha2=zeros(1,n_i);
beta2=zeros(1,n_i);
Norm=zeros(1,n_i);
for j=1:n_i
dT=0.0001;    
w2=k(j)^2+(a+b.*tanh((t-12)/dT)).^2;
w_i=sqrt(k(j)^2+c_i^2);
w_f=sqrt(k(j)^2+c_f^2);

% if rem(j,4)==0
% figure(1)
% plot(t,w2,'b')
% xlabel('t')
% ylabel('\omega (t)')
% set(gca,'fontsize',20)
% hold on
% drawnow;
% end

% RK 4

f=zeros(1,nt);
g=zeros(1,nt);

% initial condition
f(1) = sqrt(1/(2*w_i));
g(1) = -1i*(w_i)*sqrt(1/(2*w_i));

t(2) = t(1)+dt;
    f(2)=f(1)+dt*g(1)-0.5*(dt^2)*(w(1)^2)*f(1);

% new set of initial conditions is f(1) and f(2)

for n=2:nt-1
    
f(n+1) = (2*(1-(5/12)*(dt^2)*(w2(n)))*f(n)-(1+(1/12)*(dt^2)*(w2(n-1)))*f(n-1))/(1+(1/12)*(dt^2)*(w2(n+1)));

end

f(nt)=f(nt-1);
alpha=zeros(1,nt);
beta=zeros(1,nt);

for i=2:nt-1

f_out_a=sqrt(1/(2*w_f))*exp(-1i*w_f*t(i));
f_out_b=conj(f_out_a);
df_out_a=-1i*w_f*f_out_a;
df_out_b=1i*w_f*f_out_b;

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

Norm(j)=abs(alpha(nt))^2-abs(beta(nt))^2;
% alpha2(j)=abs(alpha(nt-10000))^2;
beta2(j)=abs(beta(nt))^2;

end

figure(2)
plot(k,Norm)
xlabel('k')
title('norm check ','Fontsize', 20)
set(gca,'fontsize',20)
figure(3)
plot(k,beta2)
xlabel('k')
ylabel('|\beta_{k}|^2, particle density')
title('|\beta_{k}|^2 vs k ','Fontsize', 20)
set(gca,'fontsize',20)
% figure(4)
% plot(k,alpha2)
% title('abs(alpha)^2 ','Fontsize', 20)