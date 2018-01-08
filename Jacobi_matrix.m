function Jacobi = Jacobi_matrix( u )
% Calculates Jacobi matrix
J = length(u)-1; delta = zeros(J+1);
hu = 1e-5;
for i = 1:J+1
    up = u; um = u;
    up(i) = u(i) + hu; um(i) = u(i) - hu;
    f1 = right_hand_arc_length( up );
    f2 = right_hand_arc_length( um );
    for j = 1:J+1
        delta(j,i) = f1(j) - f2(j);
    end
end
Jacobi = delta/2/hu;
end