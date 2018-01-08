function u_hat = CROS( step, u )
% Performs one step of complex one-stage Rosenbrock scheme
J = length(u)-1;
Jacobi = Jacobi_matrix( u ); E = eye(J+1);
A = E - 0.5*(1+1i)*step*Jacobi; b = right_hand_arc_length(u);
x = A\b;
u_hat = u + step*real(x);
end