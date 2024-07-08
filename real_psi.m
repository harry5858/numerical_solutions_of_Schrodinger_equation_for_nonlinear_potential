%
% % Calculate the real part of the wavefunction at time t=t+delta_t,
% t+2*delta_t etc....
% % given the value at time t. Vectorise for speed.
function [R_next]= real_psi(x, nx, R_current, I_current, dt, dx, w)

R_next= zeros(1,nx);
s=dt/(2*dx^2);

for i=2:nx-1
    
 % Calculate the real part of the wavefunction at time t=t+delta_t,
 % given the value at time t. Vectorise for speed.
 
 R_next(i)=R_current(i) - s*(I_current(i+1)-2*I_current(i)+I_current(i-1))+(0.5)*dt*(w^2)*(x(i)^2).*I_current(i);

 % Boundary conditions

 R_next(1)=0;
 R_next(nx)=0;
 
end