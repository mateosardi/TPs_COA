% clc; clear all; close all;
TamanioFuente=14;
CI=1.1; color='.r';colorc='r';

% dW1 = zeros(ho,longx+1); dW2 = zeros(dd,ho+1);
% W1 = rand(ho,longx+1); W2 = rand(dd,ho+1);
% W1a=W1-0.25; W2a=W2-0.25;


tiempo_etapa=.1;
etapas=100;
xx1=0:etapas-1;%Angulo crítico 6.1*.2 rad aprox.
sx=[0;0;CI;0];
t=xx1;

% ref = pi/2*square(2*pi/50*t);
% plot(ref);

sal=zeros(4,1);   y1=zeros(1,etapas);y2=zeros(1,etapas);
vmax=0.2;u_a=0; Mw=zeros(4,etapas-1); costo=[];costo(1)=0;consigna=[];consigna(1)=0;
color='.-k';sigma_p=.01;sigma_tita=.1;
% color='.-r';sigma_p=.01;sigma_tita=.05;
% color='.-b';sigma_p=.01;sigma_tita=.02;
for k=1:etapas-1
    
    
    VX(:,k)=sx;
    entrada = [sx(1); sx(2);sx(3);sx(4)]+1*diag([sigma_p 0 sigma_tita 0])*randn(4,1);
    %Medición ruidosa
    y1(k)=entrada(3); %angulo del brazo
    y2(k)=entrada(1); %posición del carro
    X = [ entrada; 1];     % bias
    s1 = W1a * X;
    y11=[pmntanh(s1); 1];
    s2 = W2a * y11;
    consigna(k)=s2;
    sal=mopdm2(tiempo_etapa,sx,consigna(k));
    sx=sal(:,1);
    costo_p(k)=indice_g(sal,consigna(k));
    costo(k+1)=costo(k)+indice_g(sx,consigna(k));
    u_a=consigna(k);    
end

Realizaciones=15; %Cantidad de realizaciones para el Monte Carlo.
kmax=etapas;
Kx=zeros(kmax,4);
Kv=zeros(kmax,4);
t=0:kmax- 1;u=zeros(Realizaciones,kmax);
y_sal=zeros(Realizaciones,kmax);
sQ=0.005;
F_=sQ*eye(4); %Covarianza del ruido de estado Sigma=sqrt(sQ)

for trial=1:Realizaciones %Empieza el Monte Carlo
    v=randn(4,kmax);%Señales aleatorios de media nula y varianza unidad.
    w=randn(1,kmax);
    x=[0;0;.2;0];
    for ki=1:kmax-1
        X = [ x; 1]; % bias
        s1 = W1a * X;
        y1=[pmntanh(s1); 1];
        s2 = W2a * y1;
        u(trial,ki)=s2;
%         y_sal(trial,ki)=Mat_C*x+G_*w(ki);
        
        x=mopdm2(tiempo_etapa,x,u(trial,ki))+F_*v(:,ki);
        
        p(trial,ki+1)=x(1);
        p_p(trial,ki+1)=x(2);
        alfa(trial,ki+1)=x(3);
        omega(trial,ki+1)=x(4);
    end
    u(trial,ki+1)=-Kx(1,:)*[x]-Kv(1,:)*[v(:,ki)];
end

figure(1);hold on;
subplot(3,2,1);hold on;grid on; title('Velocidad ángulo','FontSize',TamanioFuente);hold on;
plot(t,mean(omega),colorc); hold on;plot(t,mean(omega)+.5*sqrt(var(omega)),colorc);plot(t,mean(omega)- .5*sqrt(var(omega)),colorc);
subplot(3,2,2);hold on;grid on;title('Ángulo','FontSize',TamanioFuente);hold on;
plot(t,mean(alfa),colorc); hold on;plot(t,mean(alfa)+.5*sqrt(var(alfa)),colorc);plot(t,mean(alfa)- .5*sqrt(var(alfa)),colorc);
subplot(3,2,3);hold on; grid on;title('Posición carro','FontSize',TamanioFuente);hold on;
plot(t,mean(p),colorc); hold on;plot(t,mean(p)+.5*sqrt(var(p)),colorc);plot(t,mean(p)- .5*sqrt(var(p)),colorc);
subplot(3,2,4);hold on; grid on;title('Velocidad carro','FontSize',TamanioFuente);hold on;
plot(t,mean(p_p),colorc); hold on;plot(t,mean(p_p)+.5*sqrt(var(p_p)),colorc);plot(t,mean(p_p)- .5*sqrt(var(p_p)),colorc);
subplot(3,1,3); grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
plot(t,mean(u),colorc); hold on;plot(t,mean(u)+.5*sqrt(var(u)),colorc);plot(t,mean(u)- .5*sqrt(var(u)),colorc);


NPD=costo(etapas);
y1(k+1)=sal(3,1);
y2(k+1)=sx(1,1); %angulo del brazo
figure(10);
uo=consigna; To=5*6/109;
subplot(4,1,1),plot(To*xx1,y1,color),title('Evolución con neurocontrolador'),ylabel('Ángulo'),grid on,hold on;
subplot(4,1,2),plot(To*xx1,y2,color),ylabel('Posición'),grid on,hold on;
subplot(4,1,3),plot(To*xx1,costo,color),ylabel('Costo'),grid on,hold on; %,axis([0 etapas 0 6.1]);
subplot(4,1,4),plot(To*xx1(1:length(uo)),uo,color),ylabel('Acción de control'),xlabel('Tiempo [seg.]'),grid on,hold on;
