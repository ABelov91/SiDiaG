% Initial conditions: select a problem 0) autonomous power+log pole,
% 1) non-autonomous power+log pole, 2) mixed singularity, 3) S-regime of
% non-linear combustion
global problem_choice;
problem_choice = 2;
if( problem_choice == 0 || problem_choice == 1 ); initial_cond = [1;1]; end
if( problem_choice == 2 ); initial_cond = 1; end
if( problem_choice == 3 )
    t0 = 1; J = 100+1; LS = pi*sqrt(3); X = zeros(J,1);
    for j = 1:J
        X(j) = -LS/2 + LS*(j-1)/(J-1);
    end
    initial_cond = zeros(J,1);
    for j = 1:J
        initial_cond(j) = 0.5*sqrt(3)/sqrt(t0)*cos( pi*X(j)/LS );
    end
end
Time = 1; % expected time interval
epsilon_user = 1e-5; % required accuracy
[u, type, q, t0, epsilons] = SiDiaG( initial_cond, Time, epsilon_user );