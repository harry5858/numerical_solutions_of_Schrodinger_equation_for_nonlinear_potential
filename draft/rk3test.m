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

tmin=0;
tmax=10;
nt=200000;
dt1=(tmax-tmin)/nt;
t=linspace(tmin,tmax);
dt=2*dt1;
w=2;

V=0.5.*(w.*x).^2;

psi_i=((w/pi)^0.25).*exp(-(w/2).*x.^2);
psi_c=zeros(1,nx);

for j=2:nx-1
    
    k1=dt1*(-1i)*((psi_i(j+1)-2*psi_i(j)+psi_i(j-1))/(dx^2)+V(j)*psi_i(j));
    
    k2=dt1*(-1i)*(((psi_i(j+1)+k1/2)-2*(psi_i(j)+k1/2)+(psi_i(j-1)+k1/2))/(dx^2)+V(j)*(psi_i(j)+k1/2));
    
    k3=dt1*(-1i)*(((psi_i(j+1)+k2/2)-2*(psi_i(j)+k2/2)+(psi_i(j-1)+k2/2))/(dx^2)+V(j)*(psi_i(j)+k2/2));
    
    psi_c(j)=psi_i(j)+(1/6)*k1+(2/3)*k2+(1/6)*k3;

    psi_c(1)=0;
    
    psi_c(nx)=0;
    
end

psi_new=zeros(1,nx);

for n=1:nt
    
    for j=2:nx-1
    
    k1=dt1*(-1i)*((psi_c(j+1)-2*psi_c(j)+psi_c(j-1))/(dx^2)+V(j)*psi_c(j));
    
    k2=dt1*(-1i)*(((psi_c(j+1)+k1/2)-2*(psi_c(j)+k1/2)+(psi_c(j-1)+k1/2))/(dx^2)+V(j)*(psi_c(j)+k1/2));
    
    k3=dt1*(-1i)*(((psi_c(j+1)+k2/2)-2*(psi_c(j)+k2/2)+(psi_c(j-1)+k2/2))/(dx^2)+V(j)*(psi_c(j)+k2/2));
    
    psi_new(j)=psi_i(j)+(1/6)*k1+(2/3)*k2+(1/6)*k3;

    psi_new(1)=0;
    
    psi_new(nx)=0;
    
    end
    psi_c=psi_new;
    
    
        psi_imag=imag(psi_c);
    psi_real=real(psi_c);

    pd=(abs(psi_c)).^2;

%     N(a)=trapz(x,pd);

    %error_real(a)=norm(psi_real-psi_exact_r);

if rem(n,5000)==0

    figure(1)
    plot(x,pd,'k')  
    hold on
    plot(x,psi_imag,'b')
    plot(x,psi_real,'r')
    axis([-10 10 -2 2])
    hold off
    
    drawnow;
end
    
    
end







































