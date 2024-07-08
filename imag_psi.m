%
% % Calculate the imaginary part of the wavefunction at time
% t=t+delta_t/2,, t + 3*delta_t/2 etc
% % given the value at time t.
function [I_next]= imag_psi(x, nx, I_current, R_current, dt, dx, w)

I_next= zeros(1,nx);
s=dt/(2*dx^2);

for i=2:nx-1
    
 % Calculate the imaginary part of the wavefunction at time t=t+delta_t,
 % given the value at time t.
 I_next(i)=I_current(i) +s*(R_current(i+1)-2*R_current(i)+R_current(i-1))-(0.5)*dt*(w^2)*(x(i)^2).*R_current(i);
 
 % Boundary conditions

 I_next(1)=0;
 I_next(nx)=0;

end