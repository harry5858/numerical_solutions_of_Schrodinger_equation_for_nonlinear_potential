%error convergence check
clear;
% setting up dicrete coodinate system in x and t
for nx=[500 1000]
    
xmin=-10;
xmax=10;
dx=(xmax-xmin)./nx;
x=zeros(1,nx);
for j=1:nx
    x(j)=dx*(j-nx/2);
end
% x=linspace(xmin,xmax,nx);
tmin=0;
tmax=8;
if nx==500
nt=160000;
elseif nx==1000
        nt=160000;
elseif nx==2000
    nt=160000;
end
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end
% t=linspace(tmin,tmax,nt);

% ref_line=2.*(ones(1,nt));

% setting the time independent harmonic potential 0.5*(w*x)^2
w=2;
V=0.5.*(w.*x).^2;
%initial wave function
%insert overlap code
a=0;
psi_i=psiN(x,w,a);   % last input 0-4 corresponds to nth excitied state, 99 for coherent state

R_initial=real(psi_i);
I_initial=imag(psi_i);

I_next1=zeros(1,nx);
R_next1=zeros(1,nx);
% Initial run of Im(psi) and Re(psi)
for n=2:nx-1
    I_next1(n)=I_initial(n)+(dt/((dx^2)))*(R_initial(n+1)-2*R_initial(n)+R_initial(n-1))-(dt/2)*V(n)*R_initial(n);
    I_next1(1)=0;
    I_next1(nx)=0;
end
for n=2:nx-1
    R_next1(n)=R_initial(n)-(dt/((dx^2)))*(I_initial(n+1)-2*I_initial(n)+I_initial(n-1))+(dt/2)*V(n)*I_initial(n);   
    R_next1(1)=0;
    R_next1(nx)=0;
end
R_current=R_next1;
I_current=I_next1;
I_next=zeros(1,nx);
R_next=zeros(1,nx);
N=zeros(1,nt);
error_real=zeros(1,nt);
if nx==500
e1=zeros(1,nt);
elseif nx==1000
e2=zeros(1,nt);
end
for i=1:nt
    %real part of the exact ground state solution for error analysis
    psi_exact_r=psiexactr(x,w,t(i),a);
    
    [I_next]= imag_psi(x, nx, I_current, R_current, dt, dx, w);
    I_current=I_next;
    [R_next]= real_psi(x, nx, R_current, I_current, dt, dx, w);
    R_current=R_next;
    
    pd=R_current.^2+I_current.^2;
    %N(i) check whether the prob. density stays constant (stability check)
    N(i)=trapz(x,pd);
%     error_real(i)=norm(R_current-psi_exact_r);

    if nx==500
        e1(i)=R_current((nx/2)+10)-psi_exact_r((nx/2)+10);
    elseif nx==1000
        e2(i)=R_current((nx/2)+20)-psi_exact_r((nx/2)+20);
    end
        
    if rem(i, 400)==0
    figure(2)
    plot(x,R_current,'r')
    hold on
    plot(x,I_current,'b')
    plot(x,pd,'-k')
%     plot(x,psi_exact_r)
%     plot(x,V)
    hold off
    xlabel('x','Fontsize', 24)
    axis([-10 10 -2 2])
%     legend('real part of psi','imag. part of psi', 'prob. density')
    drawnow;
    end
end
% figure(3)
% plot(t,N)
% hold on
% title('normalization vs time')
% xlabel('t')
% ylabel('normalization')
% figure(4)
% hold on
% plot(t,error_real)
% axis([0 8 0 1e-2])
% title('norm of error of the real part of numerical psi vs time')
% xlabel('t')
% ylabel('L1 norm of e')

end
ref_line=2.*(ones(1,nt));
figure(4)
error_con=log2(abs(e1./e2));
plot(t,error_con)
hold on
plot(t,ref_line)
title('Convergence of error')
xlabel('t','Fontsize', 24)
set(gca,'fontsize',20)

