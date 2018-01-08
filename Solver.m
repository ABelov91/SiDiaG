function [u, q, t0, type] = Solver( step, N, B )
% Solves ODE system and analyses its singularities on a mesh with given
% step
global T; global U0;
global slope; global hypothesis_list; global stop;
if( nargin < 2 ); N = 1e+6; end
J = length(U0)+1; u = zeros(J+1,N+1);
for j = 1:J-1
    u(j,1) = U0(j);
end
u_temp = zeros(J+1,2);
H = length(hypothesis_list);
q_temp  = zeros(2,J-1,N); q  = zeros(J-1,N); 
t0_temp = zeros(2,J-1,N); t0 = zeros(J-1,N);
type = zeros(J-1,1); r = step*ones(J-1,1);
if( nargin == 3 )
    for n = 1:N
        u(:,n+1) = CROS(step, u(:,n)); 
        u_temp(:,1) = u(:,n); u_temp(:,2) = u(:,n+1);
        [q(:,n), t0(:,n), newton_root] = diagnostics(step, B, u_temp, r);
        r = newton_root;
    end
else
    n = 1; condition = 1;
    while( condition )
        u(:,n+1) = CROS( step, u(:,n) );
        u_temp(:,1) = u(:,n); u_temp(:,2) = u(:,n+1);
        for h = 1:H
            k = hypothesis_list(h); test_type = k*ones(J-1,1);
            [q_temp(k,:,n), t0_temp(k,:,n), newton_root] = ...
                diagnostics(step, test_type, u_temp, r);
            r = newton_root;
            if( n > 1 )
                for j = 1:J-1
                    q_slope  = abs( q_temp(k,j,n) -  q_temp(k,j,n-1))/step;
                    t0_slope = abs(t0_temp(k,j,n) - t0_temp(k,j,n-1))/step;
                    check = q_slope < slope && t0_slope < slope && ...
                           abs(t0_temp(k,j,n)) > u(J,2) - u(J,1);
                    if( check == 1 ); type(j) = k; end
                end
            end
        end
        for j = 1:J-1
            if( type(j) ~= 0 )
                q(j,n)=q_temp(type(j),j,n); t0(j,n)=t0_temp(type(j),j,n);
            end
        end
        t_min = T/stop;
        for j = 1:J-1
            if( abs(t0(j,n)) ~= 0 && abs(t0(j,n)) < t_min )
                t_min = abs(t0(j,n));                
            end
        end
        condition = u(J,n+1) < stop*t_min;
        n = n+1;
    end
    N = n; u_truncated  = zeros(J+1,N);
    for n = 1:N
        u_truncated(:,n)  = u(:,n);
    end
    u = u_truncated;
    q_truncated = zeros(J-1,N-1); t0_truncated = zeros(J-1,N-1);
    for n = 1:N-1
        q_truncated(:,n) = q(:,n); t0_truncated(:,n) = t0(:,n);
    end
    q = q_truncated; t0 = t0_truncated;
end
end