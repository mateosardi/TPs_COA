% Programación dinámica.
% Apredizaje Q. Para el Ejemplo 12.7.1.
% x(k+1)=x(k)+(2-2*x(k)+5/4*x(k)^2-1/4*x(k)^3)*u(k);
% con el funcional de costo J=sum((2+u(k))*exp(-x(k)));
%Autor JAP
clear,clc,close all;
Max_it=50;
TM=500; Mmax=25; color='.-k';
tic; du=Mmax; etapas=6; xmin=0; xmax=3;umin=-1; umax=1;
%%%Carga de datos
rand('state',0); equis1=3*(rand(TM,1)); tiempo=ceil((etapas-1)*rand(TM,1));
M_est = [tiempo,equis1]';
Au=(umax-umin)/(Mmax-1);
for i=1:Mmax
    uf(i)=umin+Au*(i-1);
end
CI=2;vfinal=1; Q = zeros(TM,du);
J=zeros(1,TM); sal(1)=CI;costo(1)=0;
Ya=zeros(1,TM);
for k=1:etapas-1
    entrada = [k; sal(k)];
    consigna(k) = pol_tab_mu1(entrada,M_est,Ya);
    sal(k+1)=mopdm(k,sal(k),consigna(k));
    costo(k+1)=costo(k)+indice(k,sal(k),consigna(k));
end
costo(k+1)=costo(k+1)+abs(sal(k+1)-vfinal);
evoluc(1)=costo(etapas); costo(etapas)
m=zeros(TM,du); Q11=zeros(Max_it,1);Q12=Q11;Q33=Q12;Q51=Q11;
for iterac=1:Max_it
    for iq=1:TM
        k=M_est(1,iq);
        x=M_est(2,iq);
        %Recorre todas las acciones
        for acc=1:du
            xy=mopdm(k,x,uf(acc));
            m(iq,acc)=m(iq,acc)+1;
            gama=.10*m(iq,acc)/(1.0+.10*m(iq,acc));
            if k<etapas-1
                [aux lugar]=min(abs(((k+1)-M_est(1,:))/max(M_est(1,:)))+abs((xy-M_est(2,:))/max(M_est(2,:))));
                %                 [aux lugar]=min(abs((xy-M_est(2,:))));
                Q(iq,acc)=(1-gama)*Q(iq,acc)+gama*(indice(k,x,uf(acc))+J(lugar));
            else
                Q(iq,acc)=(1-gama)*Q(iq,acc)+gama*(indice(k,x,uf(acc))+abs(xy-vfinal));
            end
            if (iq==1) && (acc== 1) %  k=M_est(1,1), x= M_est(2,1), uf(1);k =1,x =2.8504, u=-1
                Q11(iterac)=Q(iq,acc);
                f_gama(iterac)=gama;
            end
            if (iq==10) && (acc== 2) %k=M_est(1,10), x= M_est(2,10), uf(2);k =5,x =1.3341,u=-0.8571
                Q12(iterac)=Q(iq,acc);
            end
            if( iq==30) && (acc== 3) %k=M_est(1,30), x= M_est(2,30), uf(3); k =5,x =0.5964,u=-0.7143
                Q33(iterac)=Q(iq,acc);
            end
            if (iq==50) && (acc== 12) %k=M_est(1,50), x= M_est(2,50), uf(12); k =5,x =0.5690,u=0.5714
                Q51(iterac)=Q(iq,acc);
            end
        end
    end
    for iq=1:TM
        [val lugar]=min(Q(iq,:));
        J(iq)=val;
        Ya(iq)=uf(lugar);
    end
    sal(1)=CI;
    costo=0;
    for k=1:etapas-1
        entrada = [k; sal(k)];
        consigna(k) = pol_tab_mu1(entrada,M_est,Ya);
        sal(k+1)=mopdm(k,sal(k),consigna(k));
        costo(k+1)=costo(k)+indice(k,sal(k),consigna(k));
    end
    costo(k+1)=costo(k+1)+abs(sal(k+1)-vfinal);
    evoluc(iterac)=costo(etapas); costo(etapas)
    %     plot(sal),pause
end
toc
uo=consigna;
xx=0:etapas-1;
subplot(3,1,1),plot(xx,sal,color),title('Estados'),ylabel('Estados'),grid on,hold on,axis([0 5.1 0 3.1]);
subplot(3,1,2),plot(xx,costo,color),title('Costo'),ylabel('Costos'),grid on,hold on,axis([0 5.1 0 2.1]);
subplot(3,1,3),plot(xx(1:length(uo)),uo,color),title('Acción de control'),ylabel('Acciones de control'),xlabel('Etapas'),grid on,hold on,axis([0 5.1 -1.1 1.1]);
figure;
semilogy(evoluc,'.k');
xlabel('Iteraciones');
title('Costo para ir desde x(0)=2 a x(5)=1');
costo=zeros(4,etapas);estado=costo;u_opt=costo;
for CI=1:4
    in=CI-1;    estado(CI,1)=in;
    for k=1:etapas-1
        entrada = [k; in];
        an = pol_tab_mu1(entrada,M_est,Ya);
        u_opt(CI,k)=an;
        costo(CI,k+1)=indice(k,in,an)+costo(CI,k);
        in=mopdm(k,in,an);
        estado(CI,k+1)=in;
    end
    costo(CI,k+1)=costo(CI,k+1)+abs(in-1);
end
figure;
xx=0:etapas-1;
color='.-';
subplot(2,2,1),plot(xx,estado,color),title('Estados'),ylabel('Estados'),grid on;
subplot(2,2,2),plot(xx,costo,color),title('Costo'),ylabel('Costos'),grid on;
subplot(2,2,3),plot(xx(1:length(u_opt(1,:))),u_opt,color),title('Acción de control'),ylabel('Acciones de control'),xlabel('Etapas'),grid on,hold on,axis([0 5.1 -1.1 1.1]);
subplot(2,2,4),semilogy(evoluc),title('Evolucion del J(0)');xlabel('Iteraciones');
figure;
subplot(2,2,1),plot(Q11,color);
title(['Evolucion de Q_{i,u}:k=' num2str(M_est(1,1)) ', x=' num2str(M_est(2,1)) ', u=' num2str(uf(1))])
subplot(2,2,2),plot(Q12,color);
title(['Evolucion de Q_{i,u}:k=' num2str(M_est(1,10)) ', x=' num2str(M_est(2,10)) ', u=' num2str(uf(2))]);
subplot(2,2,3),plot(Q33,color);
title(['Evolucion de Q_{i,u}:k=' num2str(M_est(1,30)) ', x=' num2str(M_est(2,30)) ', u=' num2str(uf(3))]);
xlabel('Iteraciones');
subplot(2,2,4),plot(Q51,color);
title(['Evolucion de Q_{i,u}:k=' num2str(M_est(1,50)) ', x=' num2str(M_est(2,50)) ', u=' num2str(uf(12))]);
xlabel('Iteraciones');