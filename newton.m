function [root, delta, l] = newton( f0, f1, initial_approx )
% Newton method for solving scalar non-linear algebraic equation
global eps;
iter_num = 20; r = zeros(1,iter_num);
delta = 1; l = 1; r(1) = initial_approx;
while ( delta > min(1e-8,eps) )
    r(l+1) = r(l) - f0( r(l) )/f1( r(l) );
    delta = abs(r(l+1) - r(l));
    if (l > iter_num); break; end
    l = l+1;
end
root = r(l);
end