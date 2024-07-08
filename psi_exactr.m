function psi_exact_r=psiexactr(x,w,t,n)
if n==0
   psi_exact_r=((w/pi)^0.25).*exp(-(w/2).*x.^2)*cos((w/2).*t);
elseif n==1
    psi_exact_r=((w/pi)^0.25)*(sqrt(2*w)).*x.*exp(-(w/2).*x.^2).*cos((3*w/2).*t);
elseif n==2
    psi_exact_r=((w/pi)^0.25)*(1/sqrt(2)).*(2*w.*x.^2-1).*exp(-(w/2).*x.^2).*cos((5*w/2).*t);
elseif n==3
    psi_exact_r=((w/pi)^0.25)*(1/sqrt(48)).*(8*(w^1.5).*x.^3-12*sqrt(w)*x).*exp(-(w/2).*x.^2).*cos((7*w/2).*t);
elseif n==4
    psi_exact_r=((w/pi)^0.25)*(1/sqrt(384)).*(16*(w^2).*x.^4-48*(w)*x.^2+12).*exp(-(w/2).*x.^2).*cos((9*w/2).*t);
end

