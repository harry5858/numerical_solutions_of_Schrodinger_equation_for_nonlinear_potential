%evolution of Coherent state and its P and X expectation value
clear;
% setting up dicrete coodinate system in x and t
xmin=-10;
xmax=10;
nx=1000;
dx=(xmax-xmin)./nx;
x=zeros(1,nx);
for j=1:nx
    x(j)=dx*(j-nx/2);
end
% x=linspace(xmin,xmax,nx);
tmin=0;
tmax=8;
nt=160000;
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end
% t=linspace(tmin,tmax,nt);

% setting the time independent harmonic potential 0.5*(w*x)^2
w=2;
V=0.5.*(w.*x).^2;
%initial wave function
%insert overlap code
a=1;
psi_i=psiN(x,w,99);   % last input 0-4 corresponds to nth excitied state, 99 for coherent state

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
x_expected=zeros(1,nt);
p_expected=zeros(1,nt);
for i=1:nt
    %real part of the exact ground state solution for error analysis
%     psi_exact_r=psiexactr(x,w,t(i),a);

    [I_next]= imag_psi(x, nx, I_current, R_current, dt, dx, w);
    I_current=I_next;
    [R_next]= real_psi(x, nx, R_current, I_current, dt, dx, w);
    R_current=R_next;
    
    psi=R_current+1i.*I_current;
    x_element=zeros(1,nx);
    for j=1:nx-1
        x_element(j)=dx*(conj(psi(1+j))*x(1+j)*psi(1+j));
    end
    x_element2=sum(x_element);
    x_expected(i)=real((dx/2)*(conj(psi(1))*x(1)*psi(1))+x_element2+(dx/2)*(conj(psi(nx))*x(nx)*psi(nx)));
    
    p_element=zeros(1,nx);
    for j=3:nx-2
        p_element(j)=conj(psi(j))*(psi(j+1)-psi(j-1));
    end
    p_element2=sum(p_element);
    p_expected(i)=real((-1i/2)*(conj(psi(1))*(psi(2)-psi(1))+p_element2+conj(psi(nx))*(psi(nx)-psi(nx-1))));
    
    pd=R_current.^2+I_current.^2;
    
    if rem(i, 400)==0
    figure(1)
    plot(x,R_current,'r')
    hold on
    plot(x,I_current,'b')
    plot(x,pd,'-k')
    hold off
%     plot(x,psi_exact_r)
%     plot(x,V)
    xlabel('x','Fontsize', 24)
    title('real and imaginary part of psi from t=0 to t=8 for cohernet state','Fontsize', 24)
    set(gca,'fontsize',20)
    axis([-10 10 -2 2])
%     legend('real part of psi','imag. part of psi')
    figure(2)
%     plot(t,x_expected,'k')
%     hold on
%     plot(t,p_expected,'r')
%     hold off
%     grid on
%     xlabel('t','Fontsize', 24)
%     title('Momenetum expectation value and position expectation value of cohernet state from t=0 to t=8','Fontsize', 24)
%     set(gca,'fontsize',20)
%     axis([0 8 -2.1 2.1])
%     legend('position expectation value','Momenetum expectation value')
    plot(x_expected(i),p_expected(i),'k*')
    hold on
    axis([-1.5 1.5 -3 3])
    grid on
    xlabel('position expectation value','Fontsize', 24)
    ylabel('Momenetum expectation value','Fontsize', 24)
    title('Momenetum expectation value and position expectation value of cohernet state from t=0 to t=8','Fontsize', 24)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    set(gca,'fontsize',20)

%     figure(3)
%     plot(t,p_expected)
%     xlabel('t')
%     title('Momenetum expectation value of cohernet state from t=0 to t=8')
%     axis([0 8 -2.1 2.1])
    drawnow
    end
    
end
