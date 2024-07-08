function psi_i=psiNoverlap(x,w_f,n)

if n==0
   psi_i=((w_f/pi)^0.25).*exp(-(w_f/2).*x.^2);
elseif n==1
    psi_i=((w_f/pi)^0.25)*(sqrt(2*w_f)).*x.*exp(-(w_f/2).*x.^2);
elseif n==2
    psi_i=((w_f/pi)^0.25)*(1/sqrt(2)).*(2*w_f.*x.^2-1).*exp(-(w_f/2).*x.^2);
elseif n==3
    psi_i=((w_f/pi)^0.25)*(1/sqrt(48)).*(8*(w_f^1.5).*x.^3-12*sqrt(w_f)*x).*exp(-(w_f/2).*x.^2);
elseif n==4
    psi_i=((w_f/pi)^0.25)*(1/sqrt(384)).*(16*(w_f^2).*x.^4-48*(w_f)*x.^2+12).*exp(-(w_f/2).*x.^2);
end