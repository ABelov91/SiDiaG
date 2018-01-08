function f = right_hand_arc_length( u )
% Introduces arc length argument after trivial autonomization
J = length(u) - 1;
f = right_hand_autonomization( u ); f(J+1) = 1;
S = sqrt( sum(f.^2) );
f = f/S;
end