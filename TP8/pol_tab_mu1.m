%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = pol_tab_mu1(entrada,M_est,Ya);
%Look-up-table de distancia absoluta ponderada
%Da la accion de control correspondiente a la entrada
%entrada = [x k]
%Ya = [u(:)];
%M_est= [k(:) ix1(:) ix2(:)];
[a in]=min(abs((entrada(1)-M_est(1,:))/3)+abs((entrada(2)-M_est(2,:)))/6);
u=Ya(in);
%Fin rutina pol_tab_mu1
