function u = pol_tab_mu1(entrada,M_est,Ya);
%Da la accion de control correspondiente a la entrada
%entrada = [ix(1) ix(2)]
%Ya = [u(:)];
%M_est= [k(:) ix1(:) ix2(:)];
[a in]=min(abs((entrada(1)-M_est(1,:))/3)+abs((entrada(2)-M_est(2,:)))/6);
u=Ya(in);