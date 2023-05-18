% Programación dinámica.
% Apredizaje Q. Para péndulo inestable y estable
% con el funcional de costo J=sum((2+u(k))*exp(-x(k)));
%Autor JAP
clear,clc,close all;

% PENDULO EN EL EQUILIBRIO INESTABLE

% Condidiones iniciales
% x = [delta(i); delta_p(i); fi(i); fi_p(i)]
x_ini = [0;0;0.1;0];
angulo(1) = x_ini(3);
vfinal=0;

% Algoritmo aprendizaje Q
Max_it=10; % Máximas iteraciones
TM=100; % Cantidad de estados
Nmax=7; % Valores que tiene x
Mmax=9; % Valores que tiene uf
t_etapa = 0.1;
color='.-k';
tic;dx=Nmax; du=Mmax; etapas=10;umin=-15; umax=15;xmin=0; xmax=3;
%%%Carga de datos
rand('state',0);

%VALORES OBTENIDOS CON EL PLANO DE FASE
delta= 1*(randn(TM,1));
delta_p = 1*(randn(TM,1));
fi = 2*(randn(TM,1));
fi_p = 2*(randn(TM,1));

M_est = [delta, delta_p, fi, fi_p]; % Matriz de estados


CI=2; longx = 4; dd = 1; NP = dx*etapas;
ho = 3;
dW1 = zeros(ho,longx+1); dW2 = zeros(dd,ho+1);
W1 = rand(ho,longx+1); W2 = rand(dd,ho+1); W1a=W1-.5;W2a=W2-.5;
NetDef = ['HHH'
    'L--'];
maxiter = 15; stop_crit = 1e-9; lambda=1; D=0;
trparms=[maxiter stop_crit lambda D]; vfinal=1;
gama(1)=1;
tic; divi=1; xx=0:etapas-1;

Ax=(xmax-xmin)/(Nmax-1); %Delta x
Au=(umax-umin)/(Mmax-1); %Delta u


for i=1:Mmax
    uf(i)=umin+Au*(i-1);
end



for iterac=1:5
    Q = 10e3*ones(1,du);
    C=zeros(1,TM);
    J=zeros(1,TM);
    pos=0;
    for posi_x=1:TM
        x_ini = M_est(posi_x,:);
        entrada = x_ini;
        %         for k=1:etapas
        %             if k<etapas
        X = [ entrada'; 1]; % bias
        s1 = W1a * X;
        y1=[pmntanh(s1); 1];
        s2 = W2a * y1;
        u = s2;
        xy=mopdm2(t_etapa,entrada,u);
        %                 if (xy>=xmin)&(xy<=xmax)
        X = [xy; 1];
        s1 = W1 * X;
        y1=[pmntanh(s1); 1];
        s2 = W2 * y1;
        y2=s2;
        C(posi_x)=y2+indice_g(entrada',u);
        entrada = xy;
        %                 end
        %             else
        %                 m(k,posi_x)=m(k,posi_x)+1;
        %                 C(k,posi_x)=abs(x(posi_x)-vfinal);
        %             end
        
        %         end
        PHI(:,posi_x)=x_ini;
        Yo(posi_x)=C(posi_x);
    end
    [W1,W2,PI_vector,iteration,lambda]=marq(NetDef,W1,W2,PHI,Yo,trparms);
    
    for posi_x=1:TM
        %         for i=1:dx
        x_ini = M_est(posi_x,:);
        entrada = x_ini;
        for acc=1:du
            xy=mopdm2(t_etapa,x_ini,uf(acc));
            %                 if (xy>=xmin)&(xy<=xmax)
            %                     [val lugar]=min(abs(x-xy));
            entrada = xy;
            X = [ entrada; 1]; % bias
            s1 = W1 * X;
            y1=[pmntanh(s1); 1];
            s2 = W2 * y1;
            y2 = s2;
            Q(acc)=indice_g(x_ini',uf(acc))+y2;
            %                 end
        end
        [val lugar]=min(Q(:));
        J(posi_x)=val;
        Yo(posi_x)=uf(lugar);
        PHI(:,posi_x)=x_ini;
        %         end
    end
    %     pos=0; PHI=0; Yo=0;
    %     for i=1:dx
    %         for k=1:etapas-1
    %
    %             pos=pos+1;
    %             PHI(1,pos)=k;
    %             PHI(2,pos)=x(i);
    %         end
    %     end
    [W1a,W2a,PI_vector,iteration,lambda]=marq(NetDef,W1a,W2a,PHI,Yo,trparms);
%     gama(iterac+1)=100/(100+1.5*iterac);
    sal(1)=CI;
    costo=0;
    entrada = [0,0,.1,0]';
    for k=1:etapas-1
        X = [ entrada; 1]; % bias
        s1 = W1a * X;
        y1=[pmntanh(s1); 1];
        s2 = W2a * y1;
%         consigna(k)=s2;
        xy=mopdm2(t_etapa,entrada,s2);
        costo=costo+indice_g(entrada,s2);
        entrada = xy;
    end
%     costo(k+1)=costo(k+1)+abs(sal(k+1)-vfinal);
    evoluc(iterac)=costo; costo
end
break;
%--------------------


for k=1:etapas-1
    consigna(k) = pol_tab_mu1_pend(entrada,M_est,Ya);
    entrada=mopdm_pend(t_etapa,entrada,consigna(k));
    angulo(k+1) = entrada(3);
    costo(k+1)=costo(k)+indice_pend(k,entrada,consigna(k));
    %     entrada = [k; sal(k)];
end
costo(k+1)=costo(k+1);
evoluc(1)=costo(etapas); costo(etapas)
m=zeros(TM,du); Q11=zeros(Max_it,1);Q12=Q11;Q33=Q12;Q51=Q11;
for iterac=1:Max_it
    for iq=1:TM
        x=M_est(:,iq);
        %Recorre todas las acciones
        for acc=1:du
            xy=mopdm_pend(t_etapa,x,uf(acc));
            m(iq,acc)=m(iq,acc)+1;
            gama=10*m(iq,acc)/(1.0+10*m(iq,acc));
            %             if k<etapas-1
            [aux lugar]=min(  abs( xy(1) - M_est(1,:)) /max(M_est(1,:))+abs(xy(2)-M_est(2,:))/max(M_est(2,:))...
                +abs(xy(3)-M_est(3,:))/max(M_est(3,:))+abs(xy(4)-M_est(4,:))/max(M_est(4,:)));
            
            Q(iq,acc)=(1-gama)*Q(iq,acc)+gama*(indice_pend(k,x,uf(acc))+J(lugar));
            %             else
            %                 Q(iq,acc)=(1-gama)*Q(iq,acc)+gama*(indice_pend(k,x,uf(acc)));
            %             end
            if (iq==1) && (acc== 1)
                Q11(iterac)=Q(iq,acc);
                f_gama(iterac)=gama;
            end
            if (iq==10) && (acc== 2) %k=M_est(1,10), x= M_est(2,10), uf(2);k =5,x =1.3341,u=-0.8571
                Q12(iterac)=Q(iq,acc);
            end
            if( iq==30) && (acc== 3) %k=M_est(1,30), x= M_est(2,30), uf(3); k =5,x =0.5964,u=-0.7143
                Q33(iterac)=Q(iq,acc);
            end
            if (iq==50) && (acc== 11) %k=M_est(1,50), x= M_est(2,50), uf(12); k =5,x =0.5690,u=0.5714
                Q51(iterac)=Q(iq,acc);
            end
        end
    end
    for iq=1:TM
        [val lugar]=min(Q(iq,:));
        J(iq)=val;
        Ya(iq)=uf(lugar);
    end
    entrada = x_ini;
    costo=0;
    for k=1:etapas-1
        consigna(k) = pol_tab_mu1_pend(entrada,M_est,Ya);
        entrada=mopdm_pend(t_etapa,entrada,consigna(k));
        angulo(k+1) = entrada(3);
        costo(k+1)=costo(k)+indice_pend(k,entrada,consigna(k));
    end
    costo(k+1)=costo(k+1);
    evoluc(iterac)=costo(etapas); costo(etapas);
    %     plot(sal),pause
end
toc
uo=consigna;
xx=0:etapas-1;
subplot(3,1,1),plot(xx,angulo,color),title('Ángulo \phi'),ylabel('\phi'),grid on,hold on;
% axis([0 5.1 0 3.1]);
subplot(3,1,2),plot(xx,costo,color),title('Costo'),ylabel('Costos'),grid on,hold on;
% axis([0 5.1 0 2.1]);
subplot(3,1,3),plot(xx(1:length(uo)),uo,color),title('Acción de control'),ylabel('Acciones de control'),xlabel('Etapas'),grid on,hold on;
% axis([0 5.1 -1.1 1.1]);
figure;
semilogy(evoluc,'.k');
xlabel('Iteraciones');
title('Costo para ir desde x(0)=2 a x(5)=1');
costo=zeros(4,etapas);estado=costo;u_opt=costo;


figure;
subplot(2,2,1),plot(Q11,color);
title(['Evolucion de Q_{i,u}:k=' num2str(M_est(1,1)) ', x=' num2str(M_est(2,1)) ', u=' num2str(uf(1))])
subplot(2,2,2),plot(Q12,color);
title(['Evolucion de Q_{i,u}:k=' num2str(M_est(1,10)) ', x=' num2str(M_est(2,10)) ', u=' num2str(uf(2))]);
subplot(2,2,3),plot(Q33,color);
title(['Evolucion de Q_{i,u}:k=' num2str(M_est(1,30)) ', x=' num2str(M_est(2,30)) ', u=' num2str(uf(3))]);
xlabel('Iteraciones');
subplot(2,2,4),plot(Q51,color);
title(['Evolucion de Q_{i,u}:k=' num2str(M_est(1,50)) ', x=' num2str(M_est(2,50)) ', u=' num2str(uf(6))]);
xlabel('Iteraciones');
