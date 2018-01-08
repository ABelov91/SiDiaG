function [q, t0, current_root] = diagnostics(step, type, u, initial_approx)
% Calculates singularity characteristics
dim = size(u); J = dim(1) - 1;
q  = zeros(J-1,1); t0 = zeros(J-1,1); u_temp = zeros(J,2);
for j = 1:J
    u_temp(j,:) = u(j,:);
end
t = u(J+1,:); tau = t(2) - t(1);
f(:,1) = right_hand_autonomization( u(:,1) );
f(:,2) = right_hand_autonomization( u(:,2) );
current_root = step*ones(J-1,1);
for j = 1:J-1
    if( type(j) == 1 )
        q(j)  = tau/( u(j,1)/f(j,1) - u(j,2)/f(j,2) ); 
        t0(j) = q(j)*u(j,1)/f(j,1) + t(1);
    elseif( type(j) == 2 )
        a1 = f(j,1)*u(j,2); a2 = f(j,2)*u(j,1);
        f0 = @(x)( a1*(x+tau)*log(x+tau) - a2*x*log(x) );
        f1 = @(x)( a1*( log(x+tau) + 1 ) - a2*( log(x) + 1 ) );
        current_root(j) = newton( f0, f1, initial_approx(j) );
        t0(j) = t(2) + current_root(j);
        q(j) = -f(j,1)/u(j,1)*(tau + current_root(j)) ...
            *log(tau + current_root(j));
    end
end
end