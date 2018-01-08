function f = right_hand_autonomization( u )
% Performs trivial autonomization: t becomes a new unknown function
J = length(u) - 2; u_temp = zeros(J,1);
for j = 1:J
    u_temp(j) = u(j);
end
t = u(J+1);
f = ones(J+1,1); f0 = right_hand( u_temp, t );
for j = 1:J
    f(j) = f0(j);
end
f(J+1) = 1;
end