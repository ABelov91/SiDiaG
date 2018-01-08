function [u, type, q_output, t0_output, epsilon_output] =...
    SiDiaG( u0, Time, epsilon_user )
% Solves given ODE system and analyses its singularities. Solution and
% singularity characteristics are calculated with guaranteed accuracy. Mesh
% is thickened until obtained accuracy is less than user-defined error.
global T; global eps; global U0; 
global slope; global hypothesis_list; global stop;
T = Time; eps = epsilon_user; U0 = u0;
slope = 1e-5; %adjusting parameters
step  = T/100;
hypothesis_list = [1,2];
stop = 0.85;
max_grid_num = 5;
precision_u  = zeros(1,max_grid_num); 
precision_q  = zeros(1,max_grid_num);
precision_t0 = zeros(1,max_grid_num);
check = 1; m = 1;
while ( check )
    if (m > max_grid_num); break; end
    if ( m == 1)        
        [u, q, t0, type] = Solver( step );
        dim = size(u); N0 = dim(2) - 1; N = N0;
        u1 = u; q1 = q(:,N); t01 = t0(:,N);
    else
        u1 = u2; q1 = q2; t01 = t02;
    end
    m = m+1; step = step/2; N = N*2;
    [u, q, t0] = Solver( step, N, type );
    u2 = u; q2 = q(:,N); t02 = t0(:,N);
    p = 2;
    precision_pointwise_u  = richardson( u1,  u2,  p );
    precision_pointwise_q  = richardson( q1,  q2,  p );
    precision_pointwise_t0 = richardson( t01, t02, p );
    J = length(u0); u_temp = zeros(J,N+1);
    for j = 1:J
        u_temp(j,:) = u(j,:);
    end
    t = u(J+1,:);
    for n = 1:N/2+1
        f = right_hand( u_temp(:,2*n-1), t(2*n-1) );
        for j = 1:J
            precision_pointwise_u(j,n) = precision_pointwise_u(j,n)...
                                       - precision_pointwise_u(J+1,n)*f(j);
        end
    end
    precision_u(m-1)  = max(max(abs(precision_pointwise_u)));   
    precision_q(m-1)  = max(max(abs(precision_pointwise_q)));
    precision_t0(m-1) = max(max(abs(precision_pointwise_t0)));
    check = precision_u(m-1)  > epsilon_user...
         || precision_q(m-1)  > epsilon_user...
         || precision_t0(m-1) > epsilon_user;
end
M = m-1; precision_output = zeros(M,3);
for m = 1:M
    precision_output(m,1) = precision_u(m);
    precision_output(m,2) = precision_q(m);
    precision_output(m,3) = precision_t0(m);
end
epsilon_output = precision_output(M,:);
q_output  =  q(:,end); t0_output = t0(:,end);
illustrations( step, u, q, t0, precision_output );
end