function u = ley_control_tab(entrada,M_est,Ya)
% function u=ley_control(sx,W1a,W2a)
%Da la accion de control correspondiente a la entrada
%entrada = [ix(1) ix(2) ix(3) ix(4)]
%Ya = [u(:)];
%M_est= [ix1(:) ix2(:) ix3(:) ix4(:)];
alfa=entrada(3);
if abs(alfa)>pi
    while abs(alfa)>2*pi
        alfa=sign(alfa)*(abs(alfa)-2*pi);
    end
    if alfa<-pi
        alfa=2*pi+(alfa);%Lo hace positivo
    end
    if alfa>pi
        alfa=-2*pi+(alfa);%Lo hace negativo
    end
end
entrada(3)=alfa;
[a in]=min(abs((entrada(1)-M_est(1,:))/max(M_est(1,:)))...
          +abs((entrada(2)-M_est(2,:))/max(M_est(2,:)))...
          +abs((entrada(3)-M_est(3,:))/max(M_est(3,:)))...
          +abs((entrada(4)-M_est(4,:))/max(M_est(4,:))));
% [a in]=min(abs((entrada(1)-M_est(1,:)))...
%     +abs((entrada(2)-M_est(2,:)))...
%     +abs((entrada(3)-M_est(3,:)))...
%     +abs((entrada(4)-M_est(4,:))));
u=Ya(in);