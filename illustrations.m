function illustrations( step, u, q, t0, precision_output )
% Auxiliary graphics: solution and singularuty characteristics on the
% thickest mesh and accuracy curve
dim1 = size(precision_output); mesh_num = dim1(1);
dim2 = size(u); J = dim2(1)-1; N_last = dim2(2);
N0 = (N_last-1)/2^mesh_num; N = zeros(1,mesh_num);
for m = 1:mesh_num
    N(m) = N0*2^(m);
end
u_temp = zeros(J-1,N_last);
for j = 1:J-1
    u_temp(j,:) = u(j,:);
end
t = u(J,:);
figure; hold on; xlabel('t'); ylabel('u'); title('SOLUTION');
for j = 1:J-1
    c = j/5 - floor(j/5);
    plot(t,u_temp(j,:),'Color',[c,c,c],'LineWidth',2);
end
str = 'u(1)';
for j = 1:J-2
    str0 = strcat('u(',mat2str(j+1),')'); str = char(str,str0);
end
if( J<=10 ); legend(str); end
figure; hold on; xlabel('l'); ylabel('q, t0'); 
title('SINGULARITY CHARACTERISTICS');
l = 0:step:(N_last-2)*step;
for j = 1:J-1
    c = j/5 - floor(j/5);
    plot(l, q(j,:),'-','Color',[c,c,c],'LineWidth',2);
    plot(l,t0(j,:),':','Color',[c,c,c],'LineWidth',2.5);
end
str = char( 'q(1)', 't0(1)' );
for j = 1:J-2
    str_q  = strcat('q(', mat2str(j+1),')');
    str_t0 = strcat('t0(',mat2str(j+1),')');
    str = char( str, str_q ); str = char( str, str_t0 );
end
if( J<=10 ); legend(str); else legend('q','t0'); end
figure; hold on; xlabel('lg N'); ylabel('lg epsilon'); title('CONVERGENCE')
plot(log10(N),log10(precision_output(:,1)),...
    '-ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',14)
plot(log10(N),log10(precision_output(:,2)),...
    '-^k','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',14)
plot(log10(N),log10(precision_output(:,3)),...
    '-sk','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',14)
legend('u','q','t0');
end