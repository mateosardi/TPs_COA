clear all; clc; close all;
randn('state',100);
dbclear if infnan
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
TamanioFuente=14;
%Condiciones iniciales
alfa(1)=.1; color='.r';colorc='r';
% alfa(1)=.5; color='.g';colorc='g';
% alfa(1)=.8; color='.b';colorc='b';
Realizaciones=15; %Cantidad de realizaciones para el Monte Carlo.
%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_Ac=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0];
Mat_Bc=[0; 1/M; 0; -1/(long*M)];
kmax=2000;
sQ=0.005;
F_=sQ*eye(4); %Covarianza del ruido de estado Sigma=sqrt(sQ)
sR=.01; %Covarianza del ruido de medicion sigma=sqrt(sR)
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
Q=diag([1e1 1e1 1e1 1e1]);%Matrices de diseño del controlador DLQG
S=diag([1e1 1e3 1e6 1e2]);
P=S; %condición inicial de P
R=1e0;
Kx=zeros(kmax,4);
Kv=zeros(kmax,4);
for hi=kmax-1:-1:1
P= Q + Aa'*P*Aa - Aa'*P*Ba*inv(R+Ba'*P*Ba)*Ba'*P*Aa;
Kx(hi,:)=inv(R+Ba'*P*Ba)*Ba'*P*Aa;
Kv(hi,:)=inv(R+Ba'*P*Ba)*Ba'*P*F_;
Ea(:,hi)=eig(Aa-Ba*Kx(hi,:));
end
Jmin=x0'*P*x0;J=0;t=0:kmax- 1;u=zeros(Realizaciones,kmax);Jn_=zeros(Realizaciones,kmax);
y_sal=zeros(Realizaciones,kmax);
for trial=1:Realizaciones %Empieza el Monte Carlo
v=randn(4,kmax);%Señales aleatorios de media nula y varianza unidad.
w=randn(1,kmax);
x=x0;
for ki=1:kmax-1
u(trial,ki)=-Kx(ki,:)*x-Kv(ki,:)*v(:,ki);
Jn_(trial,ki+1)=Jn_(trial,ki)+(x'*Q*x + u(trial,ki)'*R*u(trial,ki));
y_sal(trial,ki)=Mat_C*x+G_*w(ki);

x=mopdm2(Ts,x,u(trial,ki))+F_*v(:,ki);

p(trial,ki+1)=x(1);
p_p(trial,ki+1)=x(2);
alfa(trial,ki+1)=x(3);
omega(trial,ki+1)=x(4);
end
Jn_(trial,ki+1)=Jn_(trial,ki+1)+x'*S*x;
u(trial,ki+1)=-Kx(1,:)*[x]-Kv(1,:)*[v(:,ki)];
end

% Gráficos
t=t*Ts;Jn=mean(Jn_);disp(['El valor de costo es Jn(end)=' num2str(Jn(end)) '. Alfa(1)=' num2str(alfa(1)) '.']);
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
figure(2);hold on;
subplot(2,2,1);hold on; plot(alfa,omega,color); grid on;xlabel('Ángulo','FontSize',TamanioFuente);ylabel('Velocidad angular','FontSize',TamanioFuente);hold on;
subplot(2,2,2);hold on; plot(p,p_p,color); grid on;xlabel('Posición carro','FontSize',TamanioFuente);ylabel('Velocidad carro','FontSize',TamanioFuente);hold on;
subplot(2,2,3);hold on;
plot(t,Jn,color);plot(t,Jmin*ones(size(t)),colorc);ylabel('Acumulación de costo','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);
subplot(2,2,4);hold on; plot(abs(Ea)');ylabel('Polos de lazo cerrado','FontSize',TamanioFuente);xlabel('Etapas de iteración','FontSize',TamanioFuente);
