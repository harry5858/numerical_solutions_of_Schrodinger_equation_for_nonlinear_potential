clear;
clear;
xmin=-10;
xmax=10;
nx=1000;
dx=(xmax-xmin)./nx;
% x=zeros(1,nx);
% for j=1:nx
%     x(j)=dx*(j-nx/2);
% end
x=linspace(xmin,xmax,nx);
tmin=0;
tmax=12;
nt=240000;
dt=(tmax-tmin)/nt;
% t=zeros(1,nt);
% for i=1:nt
%     t(i)=i*dt;
% end
t=linspace(tmin,tmax,nt);
dT=zeros(1,17);
for j=1:17
 
   w_i=3;
   w_f=6;
   alpha=(w_i+w_f)/2;
   beta=(w_f-w_i)/2;
   T=tmax/2;
   dT=0.2:0.1:1.8;
   w=alpha+beta.*(tanh((t-T)./dT(j)));
   figure(1)
   plot(t,w,'b')
   hold on
   xlabel('t','Fontsize', 24)
   ylabel('\omega','Fontsize', 24)
   title('\omega verus time with different \Delta T','Fontsize', 24)
   axis([0 tmax 2 7])


   V_i=0.5*(w_i.*x).^2;
   
   psi_i=((w_i/pi)^0.25).*exp(-(w_i/2).*x.^2);
   
   R_initial=real(psi_i);
   I_initial=imag(psi_i);

   I_next1=zeros(1,nx);
   R_next1=zeros(1,nx);
   % Initial run of Im(psi) and Re(psi)
   for n=2:nx-1
    I_next1(n)=I_initial(n)+(dt/(2*(dx^2)))*(R_initial(n+1)-2*R_initial(n)+R_initial(n-1))-(dt/2)*V_i(n)*R_initial(n);
    I_next1(1)=0;
    I_next1(nx)=0;
   end
   for n=2:nx-1
    R_next1(n)=R_initial(n)-(dt/(2*(dx^2)))*(I_initial(n+1)-2*I_initial(n)+I_initial(n-1))+(dt/2)*V_i(n)*I_initial(n);   
    R_next1(1)=0;
    R_next1(nx)=0;
   end
   
   R_current=R_next1;
   I_current=I_next1;
   I_next=zeros(1,nx);
   R_next=zeros(1,nx);
   
   for i=1:nt
    
    V=0.5*(w(i).*x).^2;
    
    [I_next]= imag_psi(x, nx, I_current, R_current, dt, dx, w(i));
    I_current=I_next;
    [R_next]= real_psi(x, nx, R_current, I_current, dt, dx, w(i));
    R_current=R_next;

    pd=R_current.^2+I_current.^2;
    
    if rem(i, 5000)==0
        
    figure(2)
    plot(x,R_current,'r')
    hold on
    plot(x,I_current,'b')
    plot(x,pd,'-k')
    plot(x,V)
    hold off
    xlabel('x','Fontsize', 24)
    title('Probability density, real and imaginary part of psi','Fontsize', 24)
    set(gca,'fontsize',20)
    axis([-10 10 -2 2])
    drawnow;
    
    end
    
    end
   psi_random_1=R_current+1i.*I_current;

   pd_random_1=(abs(psi_random_1)).^2;
   %A is the normalization constant
   A=1./sqrt((trapz(x,pd_random_1)));
   %normalized psi
   psi_random_2=A.*psi_random_1;
   %normalization check
   pd_random_2=(abs(psi_random_2)).^2;
   B=trapz(x,pd_random_2); % B should be 1 if psi_random is normalized

   %c(i) is square root of the probability amplitude in the i-th eigenstate
   c=zeros(1,5);
   %C(i) is the probability amplitude in the i-th eigenstate
   C=zeros(1,5);
   for a=1:5
    
    overlap=conj(psiNoverlap(x,w_f,a-1)).*psi_random_2;
    
    c(a)=trapz(x,overlap);
    C(a)=(abs(c(a)))^2;
    
%     disp(C(i))
   end
   figure(4)
   a1=C(1);
   a2=C(2);
   a3=C(3);
   a4=C(4);
   a5=C(5);
   semilogy(dT(j),a1,'b*','MarkerSize', 15)
   hold on
   semilogy(dT(j),a2,'k*','MarkerSize', 15)
   semilogy(dT(j),a3,'ro','MarkerSize', 15)
   semilogy(dT(j),a4,'g+','MarkerSize', 15)
   semilogy(dT(j),a5,'md','MarkerSize', 15)
   axis([0 1.6 0 1])
   xlabel('\Delta T','Fontsize', 24)
   title('probabilities of the final wave function in terms of final state eigenfunctions in different \Delta T','Fontsize', 24)
   legend('C0','C1','C2','C3','C4')
   set(gca,'fontsize',20)
   
end
