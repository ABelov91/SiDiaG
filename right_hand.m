function f = right_hand ( u,t )
global problem_choice;
if( problem_choice == 0 ) % autonomous problem
    q = 2;
    f(1) = (q*( u(1)^(1-1/q))*exp(u(1)^(1/q))); % logarithmic pole
    f(2) = u(2)^3; % power pole
    f = f';
elseif( problem_choice == 1 ) % non-autonomous problem       
    q = 2.5; t0 = 0.3; % power pole
    f(1) = ( -sin(6*pi*t)/( 2 + cos(6*pi*t) ) + q/(t0 - t) )*u(1);
    f(2) = -sin(6*pi*t)/(2 + cos(6*pi*t))*u(2) ...
        + q*( u(2)^(1-1/q) )*exp( u(2)^(1/q) ); % logarithmic pole 
    f = f';
elseif( problem_choice == 2 ) % mixed singularuty
    t0 = 0.8;  q = 0.1;
    f(1) = ( -u(1)/log(t0 - t) )^( (q+1)/q ) + q*u/(t0-t);
elseif( problem_choice == 3 ) % S-regime of combustion
    J = 100+1; hx = pi*sqrt(3)/(J-1);
    f = zeros(J,1); f(1) = u(1)^3; f(J) = u(J)^3;
    for j = 2:J-1
        f(j) = 0.5*( u(j+1)^2 + u(j)^2 )*( u(j+1) - u(j) )/hx^2 ...
             - 0.5*( u(j)^2 + u(j-1)^2 )*( u(j) - u(j-1) )/hx^2 + u(j)^3;
    end
end
end