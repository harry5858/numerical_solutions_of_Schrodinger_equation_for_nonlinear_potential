function psi_i=psiN(x,w,n)

if n==0
   psi_i=((w/pi)^0.25).*exp(-(w/2).*x.^2);
elseif n==1
    psi_i=((w/pi)^0.25)*(sqrt(2*w)).*x.*exp(-(w/2).*x.^2);
elseif n==2
    psi_i=((w/pi)^0.25)*(1/sqrt(2)).*(2*w.*x.^2-1).*exp(-(w/2).*x.^2);
elseif n==3
    psi_i=((w/pi)^0.25)*(1/sqrt(48)).*(8*(w^1.5).*x.^3-12*sqrt(w)*x).*exp(-(w/2).*x.^2);
elseif n==4
    psi_i=((w/pi)^0.25)*(1/sqrt(384)).*(16*(w^2).*x.^4-48*(w)*x.^2+12).*exp(-(w/2).*x.^2);
elseif n==99
    k0=0;
    x0=1;
    psi_i=((w/pi)^0.25).*exp(-(w/2).*(x-x0).^2).*exp(-1i*k0.*x);
end
    