clear all;
for nt=[200000 100000]
tmin=0;
tmax=10;
% nt=200000;
nt1=nt/2;
t=linspace(tmin,tmax,nt);
t1=linspace(tmin,tmax,nt1);
dt=2*tmax/nt;
dt2=tmax/nt;

w_i=2;
% w_f=6;
% a=(w_f+w_i)/2;
% b=(w_f-w_i)/2;
% dT=0.0001;
% w=a+(b*tanh((t-(tmax/2))/dT));
w=2*ones(1,nt);
% RK 4
y_d=zeros(1,nt);
y1_d=zeros(1,nt);
% initial condition
y_d(1)=1/sqrt(2*w_i);
y1_d(1)=-1i*w_i*1/sqrt(2*w_i);

% k0= (dt^2)*(-1*(w(1)^2)*y(1));
% k1= (dt^2)*(-1*(w(2)^2)*(y(1)+dt2*y1(1)+k0/8));
% k2= (dt^2)*(-1*(w(3)^2)*(y(1)+dt*y1(1)+0.5*k1));
% 
% y(3)=y(1)+dt*y1(1)+(k0+2*k1)/6;
% y1(3)=y1(1)+(k0+4*k1+k2)/(6*dt);

for n=1:2:nt-2
k0= (dt^2)*(-1*(w(n)^2)*y_d(n));
k1= (dt^2)*(-1*(w(n+1)^2)*(y_d(n)+dt2*y1_d(n)+k0/8));
k2= (dt^2)*(-1*(w(n+2)^2)*(y_d(n)+dt*y1_d(n)+0.5*k1));

y_d(n+2)=y_d(n)+dt*y1_d(n)+(k0+2*k1)/6;
y1_d(n+2)=y1_d(n)+(k0+4*k1+k2)/(6*dt);

end

y=zeros(1,nt1);
y(1)=y_d(1);
for i=2:nt1
    y(i)=y_d(2*i-1);
end
figure(1)
plot(t1,real(y))
hold on
plot(t1,imag(y))
f_ex=0.5*exp(-1i*2*t);
    if nt==200000
    e=zeros(1,nt1);
    for i=1:nt1-1
    e(i)=y(i+1)-f_ex(i);
    end
%     figure(2)
%     plot(t,real(e))
%     hold on
%     drawnow;
    elseif nt==100000
        e2=zeros(1,nt1);
    for i=1:nt1-1
    e2(i)=y(i+1)-f_ex(i);
    end
    end
end
t2=linspace(tmin,tmax,100000/2);
p=zeros(1,50000);
for i=1:50000
    p(i)=log2(abs(e2(i)/e(i+1)));
end
figure(3)
plot(t1,p)

y_ex=(1/2)*exp(-1i*2*t1);
e1=y-y_ex;
figure(2)
plot(t1,e1)

% if nt==200000
%     y1=y;
% elseif nt==100000
%     y2=y;
% elseif nt==50000
%     y4=y;
% end
% end
% t2=linspace(tmin,tmax,25000);
% p_num1=zeros(1,25000);
% p_den1=zeros(1,25000);
% p1=zeros(1,25000);
% for j=1:25000
% p_num1(j) =real( y4(j)- y2(2*j));
% p_den1(j) =real( y2(2*j)- y1(4*j));
% p1(j)=p_num1/p_den1;
% end
% figure(1)
% plot(t1,log2(p1))
% alpha=zeros(1,nt1);
% beta=zeros(1,nt1);
% for i=2:nt1-1
% 
% f_out_a=sqrt(1/(2*w_f))*exp(-1i*w_f*t1(i));
% f_out_b=conj(f_out_a);
% df_out_a=-1i*w_f*f_out_a;
% df_out_b=1i*w_f*f_out_b;
% 
% df_in= (conj(y(i+1))-conj(y(i-1)))/(2*dt);
% 
% alpha(i)=1i*(conj(y(i))*df_out_a-df_in*f_out_a);
% beta(i)=1i*((conj(y(i))*df_out_b-df_in*f_out_b));
% 
% alpha(nt1-1)=alpha(nt1-2);
% alpha(nt1)=alpha(nt1-1);
% alpha(1)=alpha(2);
% beta(nt1-1)=beta(nt1-2);
% beta(nt1)=beta(nt1-1);
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
% g=zeros(1,nt1);
% 
% for j=1:nt1
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
% alpha_p=alpha_1(nt1-3);
% beta_p=beta_1(nt1-3);
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
% plot(n(1),C0,'bo','Markersize', 15)
% hold on
% plot(n(3),C2,'bo','Markersize', 15)
% plot(n(5),C4,'bo','Markersize', 15)
% txt0 = '\leftarrow |C_0|^2 = 0.9428';
% txt2 = '\leftarrow |C_2|^2 = 0.05238';
% txt4 = '\downarrow |C_4|^2 = 0.004364';
% t0=text(0,C0,txt0);
% t2=text(2,C2,txt2);
% t4=text(4,C4+0.03,txt4);
% t0.FontSize = 20;
% t2.FontSize = 20;
% t4.FontSize = 20;
% xlabel('n','Fontsize', 24)
% title('Transition probabilities |C(n)|^2 for n=0,1,2,3,4 ')
% axis([-0.1 5 0 1])
% set(gca,'fontsize',20)


