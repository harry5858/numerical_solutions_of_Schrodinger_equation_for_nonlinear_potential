clear;

tmin=0;
tmax=12;
nt=120000;
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end

c_i=3;
c_f=6;
a=(c_f+c_i)/2;
b=(c_f-c_i)/2;
dT=1;
c=a+b*tanh((t-6)/dT);

k=0;

w=sqrt(k^2+c.^2);
w_in=sqrt(k^2+c_i^2);
w_out=sqrt(k^2+c_f^2);
figure(1)
plot(t,w)

% setup for solutions
f=zeros(1,nt);
g=zeros(1,nt);

%initial condition
f(1)=sqrt(1/(2*w_in));
g(1)=-1i*w_in*sqrt(1/(2*w_in));

%initial run

f(2)=f(1)-2*dt*((w(1)^2)*g(1));
g(2)=g(1)+2*dt*f(1);

for i=2:nt-1
    
    f(i+1)=f(i-1)-2*dt*((w(i)^2)*g(i));
    g(i+1)=g(i-1)+2*dt*f(i);

    g(nt)=g(nt-1);
    f(nt)=f(nt-1);
    
end

figure(2)
plot(t,real(f))
hold on
plot(t,imag(f))

alpha=zeros(1,nt);
beta=zeros(1,nt);

for j=1:nt-1
    f_out=sqrt(1/(2*w_out))*exp(-1i*w_out*t(j));
    f_out_star=conj(f_out);
    df_out=-1i*w_out*f_out;
    df_out_star=1i*w_out*f_out_star;
    
    dg=(g(j+1)-g(j))/dt;
    
    alpha(i)=1i*(conj(g(i))*df_out-dg*f_out);
    beta(i)=-1i*(g(i)*df_out_star-dg*f_out_star);
    
    alpha(nt)=alpha(nt-1);
    beta(nt)=beta(nt-1);
    
end
% 
% figure(3)
% plot(t,abs(beta))
%     
% figure(4)
% plot(t,abs(alpha))













































































