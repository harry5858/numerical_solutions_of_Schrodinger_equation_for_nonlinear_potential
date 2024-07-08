%use ode45 to solve the equation of motion
trange=linspace(0,12,60000);
dt=max(trange)/60000;
nt=numel(trange);
w_i=3;
w_f=1:1:20;
n_i=numel(w_f);

dw=zeros(1,n_i);
for n=1:n_i
dw(n)=w_f(n)-w_i;
end

C0=zeros(1,n_i);
C2=zeros(1,n_i);
C4=zeros(1,n_i);
for j=1:n_i
dT=1/2;
[t, state_values] = ode45(@(t,x) myfunc(t,x,w_f(j),w_i,dT),trange,[sqrt(1/(2*w_i)),-1i*(w_i)*sqrt(1/(2*w_i))]);
f = state_values(:,1);
fdot = state_values(:,2);
% figure(1)
% plot(t,real(f))
% hold on
% plot(t,imag(f))
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

% Proability
% Only consider the first 3 even terms for comparison with the Schrodinger
% numerical solution
% 
alpha_p=alpha_1(nt-6000);
beta_p=beta_1(nt-6000);

A=conj(beta_p)/(2*alpha_p);
N=sqrt(1/abs(alpha_p));

C0(j)=(abs(N))^2;
C2(j)=(abs(N*sqrt(2)*(A)))^2;
C4(j)=(abs(N*sqrt(6)*(A^2)))^2;
end
figure(8)
semilogy(dw,C0,'b*','Markersize', 15)
hold on
semilogy(dw,C2,'r')
semilogy(dw,C4,'k')
% axis([-3 7.5 0 1])
legend('C0','C2','C4')
% legend(' n=2 ',' n=4 ')
xlabel('dw','Fontsize', 24)
ylabel('probability |Cn|^2','Fontsize', 24)
title('probability |<2n out | 0 in>|^2 , n=0,1,2 ','Fontsize', 20)
