clear all; clc; close all;
randn('state',1);
dbclear if infnan
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
TamanioFuente=14;
%Condiciones iniciales
alfa(1)=.1; color='.r';colorc='r';
% alfa(1)=.15; color='.g';colorc='g';
% alfa(1)=.2; color='.b';colorc='b';
Realizaciones=30; %Cantidad de realizaciones para el Monte Carlo.
%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_Ac=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0];
Mat_Bc=[0; 1/M; 0; -1/(long*M)];
kmax=2000;
sQ=1e-5;
F_=sQ*eye(4); %Covarianza del ruido de estado Sigma=sqrt(sQ)
sR=0.01; %Covarianza del ruido de medicion sigma=sqrt(sR)
G_=sR;
I=eye(4);
Mat_C=[1 0 0 0];
sys_c=ss(Mat_Ac,Mat_Bc,Mat_C,0);Ts=0.01;
sys_d=c2d(sys_c,Ts,'zoh');
Mat_A=sys_d.a;
Mat_B=sys_d.b;
Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B ];%Matriz Controlabilidad
rango=rank(Mat_M);
x=[0;0;-alfa(1);0];x0=x;
p(1)=x(1); p_p(1)=x(2); alfa(1)=x(3); omega(1)=x(4);
Aa=Mat_A;Ba=Mat_B;
%_____________ESTIMADOR KALMAN______________
P_Kalman=F_*F_';
P11=zeros(1,5000);P22=P11;P33=P11;P44=P11;
for h_k=1:5000
    P_Kalman_=Mat_A*P_Kalman*Mat_A'+(F_*F_');
    K_Kalman=P_Kalman_*Mat_C'/(Mat_C*P_Kalman_*Mat_C'+(G_*G_')); %Ganancia de Kalman
    P_Kalman=(eye(4)-K_Kalman*Mat_C)*P_Kalman_;
    P11(h_k)=P_Kalman(1,1);
    P22(h_k)=P_Kalman(2,2);
    P33(h_k)=P_Kalman(3,3);
    P44(h_k)=P_Kalman(4,4);
end
figure(4);semilogy(P11);hold on;semilogy(P22,'g');semilogy(P33,'c');semilogy(P44,'r');title('Evolución de P_1_1,P_2_2,P_3_3 y P_4_4.'); xlabel('Iteraciones');
EK=abs(eig(Mat_A-K_Kalman*Mat_C));
% break
Q=1e1*diag([1e2 1e1 1e2 1e1]);%Matrices de diseño del controlador DLQG
S=Q;
P=S; %condición inicial de P
R=1e0;
Kk=zeros(kmax,4);
for hi=kmax-1:-1:1
    P= Q + Aa'*P*Aa - Aa'*P*Ba/(R+Ba'*P*Ba)*Ba'*P*Aa;
    Kk(hi,:)=(R+Ba'*P*Ba)\Ba'*P;
    Ea(:,hi)=eig(Aa-Ba*Kk(hi,:)*Aa);
end
% plot(abs(Ea)') break
H=inv([eye(4) Ba*inv(R)*Ba'; zeros(4) Aa'])*[Aa zeros(4);-Q eye(4)];
[V,D]=eig(H);MX1X2=[];
for ii=1:8
    if abs(D(ii,ii))<1
        MX1X2=[MX1X2 V(:,ii)];
    end
end
MX1=MX1X2(1:4,:); MX2=MX1X2(5:8,:);
Pa=real(MX2*inv(MX1));% [K1,P,E]=dlqr(Mat_A,Mat_B,Q,R);
Klqr=(R+Ba'*Pa*Ba)\Ba'*Pa*Aa; % no es lo mismo que Kx(1,:)
Elqr=abs(eig(Aa-Ba*Klqr));%break
%Cálculo del observador
Qo=1e1*diag([1e1 1e0 1e1 1e0]); Ro=1e1;
Aad=Mat_A';
Bad=Mat_C';
%Contrucción del Hamiltoniano para el cálculo del controlador
Hm=([eye(4) Bad*inv(Ro)*Bad'; zeros(4) Aad'])\[Aad zeros(4);-Qo eye(4)];
[V,D]=eig(Hm);MX1X2=[];
for ii=1:8
    if abs(D(ii,ii))<1
        MX1X2=[MX1X2 V(:,ii)];
    end
end
MX1=MX1X2(1:4,:); MX2=MX1X2(5:8,:);
Po=real(MX2*inv(MX1));
Tx=(inv(Ro+Bad'*Po*Bad)*Bad'*Po*Aad)';
Em=eig(Mat_A-Tx*Mat_C); % break
Jmin=x0'*P*x0;J=0;t=0:kmax- 1;u=zeros(Realizaciones,kmax);Jn_=zeros(Realizaciones,kmax);
y_sal=zeros(Realizaciones,kmax);
for trial=1:Realizaciones %Empieza el Monte Carlo
    v=randn(4,kmax);%Señales aleatorios de media nula y varianza unidad.
    w=randn(1,kmax);p(trial,1)=x0(1);
    p_p(trial,1)=x0(2);alfa(trial,1)=x0(3);omega(trial,1)=x0(4);
    x=x0;ua=0;x_hat=zeros(size(x0));x_hat_=x_hat;
    for ki=1:kmax-1
        y_sal(trial,ki)=Mat_C*x+G_*w(ki);
        x_hat=x_hat_+K_Kalman*(y_sal(trial,ki) -Mat_C*x_hat_);
        u(trial,ki)=-Klqr*x_hat; Ley='LQR sin F';
%         u(trial,ki)=-Kk(1,:)*x_hat; Ley='LQR VT sin A ni F';
        x=mopdm2(Ts,x,u(trial,ki))+F_*v(:,ki);
        y_sal_O(trial,ki)=Mat_C*x_hat;
        Jn_(trial,ki+1)=Jn_(trial,ki)+(x'*Q*x + u(trial,ki)'*R*u(trial,ki));
        x_hat_=Mat_A*x_hat+Mat_B*u(trial,ki);
        p(trial,ki+1)=x(1);
        p_p(trial,ki+1)=x(2);alfa(trial,ki+1)=x(3);omega(trial,ki+1)=x(4);
    end
    Jn_(trial,ki+1)=Jn_(trial,ki+1)+x'*S*x;
end
t=t*Ts;Jn=mean(Jn_);disp(['Con ' Ley ' se tiene Jn(end)=' num2str(Jn(end)) '.Alfa(1)=' num2str(alfa(1)) '.']);
figure(1);hold on;
subplot(3,2,1);hold on;grid on; title('Velocidad ángulo','FontSize',TamanioFuente);hold on;
plot(t,mean(omega),colorc); hold on;plot(t,mean(omega)+.5*sqrt(var(omega)),colorc);plot(t,mean(omega)- .5*sqrt(var(omega)),colorc);
subplot(3,2,2);hold on;grid on;title('Ángulo','FontSize',TamanioFuente);hold on;
plot(t,mean(alfa),colorc); hold on;plot(t,mean(alfa)+.5*sqrt(var(alfa)),colorc);plot(t,mean(alfa)- .5*sqrt(var(alfa)),colorc);
subplot(3,2,3);hold on; grid on;title('Posición carro','FontSize',TamanioFuente);hold on;
plot(t,mean(p),colorc); hold on;plot(t,mean(p)+.5*sqrt(var(p)),colorc);plot(t,mean(p)-.5*sqrt(var(p)),colorc);
subplot(3,2,4);hold on; grid on;title('Velocidad carro','FontSize',TamanioFuente);hold on;
plot(t,mean(p_p),colorc); hold on;plot(t,mean(p_p)+.5*sqrt(var(p_p)),colorc);plot(t,mean(p_p)- .5*sqrt(var(p_p)),colorc);
subplot(3,1,3); grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
plot(t,mean(u),colorc); hold on;plot(t,mean(u)+.5*sqrt(var(u)),colorc);plot(t,mean(u)-.5*sqrt(var(u)),colorc);
figure(2);hold on;
subplot(2,2,1);hold on; plot(alfa,omega,color); grid on;xlabel('Ángulo','FontSize',TamanioFuente);ylabel('Velocidad angular','FontSize',TamanioFuente);hold on;
subplot(2,2,2);hold on; plot(p,p_p,color); grid on;xlabel('Posición carro','FontSize',TamanioFuente);ylabel('Velocidad carro','FontSize',TamanioFuente);hold on;
subplot(2,2,3);hold on;
plot(t,Jn,color);plot(t,Jmin*ones(size(t)),colorc);ylabel('Acumulación de costo','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);
subplot(2,2,4);hold on; plot(abs(Ea)');ylabel('Polos de lazo cerrado','FontSize',TamanioFuente);xlabel('Etapas de iteración','FontSize',TamanioFuente);
