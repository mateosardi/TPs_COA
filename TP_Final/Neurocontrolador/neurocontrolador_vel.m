% Programación dinámica.
%Autor JAP
clear,clc,close all;
% MOTOR DE CC

% Condidiones iniciales
% x = [ia(i); if(i) w(i); theta(i)]
x_ini = [0;0;0;0];
angulo(1) = x_ini(4);
torque = 0.1; % Torque aplicado
Vf = 220; % Tensión de campo


% Algoritmo aprendizaje Q
Max_it=15; % Máximas iteraciones
TM=1000; % Cantidad de estados
Mmax=43; % Valores que tiene uf
t_etapa = 1e-2;
color='.-k';
tic; du=Mmax;
etapas=201;
umin=-100; umax=110;
%%%Carga de datos
rand('state',0);
randn('state',0);

%VALORES OBTENIDOS CON EL CONTROLADOR PID
% Valores máximos
ia_m = 15;
if_m = 2.5;
w_m = 20;
theta_m = 2;
V_max = [1/ia_m; 1/if_m; 1/w_m; 1/theta_m]; % Vector de valores máximos
% V_max = [1;1; 1; 1];

ia= 2*ia_m*((rand(TM,1)-.5));
i_f= if_m*(ones(TM,1));
w = 2*w_m*((rand(TM,1)-.5));
theta = 2*theta_m*((rand(TM,1)-.5));

M_est = [ia,i_f, w, theta]; % Matriz de estados

CI=x_ini(4); longx = 4; dd = 1;
ho = 3;
dW1 = zeros(ho,longx+1); dW2 = zeros(dd,ho+1);
W1 = rand(ho,longx+1); W2 = rand(dd,ho+1); W1a=W1-.5;W2a=W2-.5;
NetDef = ['HHH'
    'L--'];
maxiter = 10; stop_crit = 1e-9; lambda=1; D=0;
trparms=[maxiter stop_crit lambda D];
gama(1)=1;
tic; divi=1; xx=0:etapas-1;

Au=(umax-umin)/(Mmax-1); %Delta u

for i=1:Mmax
    uf(i)=umin+Au*(i-1);
end

for iterac=1:Max_it
    Q = 10e3*ones(1,du);
    C=zeros(1,TM);
    J=zeros(1,TM);
    pos=0;
    for posi_x=1:TM
        entrada = M_est(posi_x,:)';
        PHI(:,posi_x)=entrada .* V_max;
        C_C=0;
        for k_k=1:etapas
            X = [ entrada.* V_max; 1]; % bias
            s1 = W1a * X;
            y1=[pmntanh(s1); 1];
            s2 = W2a * y1;
            u = s2;
            xy=mopdm2_motor(t_etapa,entrada,torque,u,Vf);
            xy = xy(1:4);
            C_C=C_C+indice_g_motor(entrada,u);
            entrada = xy;
        end
        Yo(posi_x)=C_C;
    end
    [W1,W2,PI_vector,iteration,lambda]=marq(NetDef,W1,W2,PHI,Yo,trparms);
    C_Yo = Yo;
    Ya = zeros(TM,1);
    for posi_x=1:TM
        x_ini = M_est(posi_x,:);
        entrada = x_ini;
        for acc=1:du
            xy=mopdm2_motor(t_etapa,x_ini,torque,uf(acc),Vf);
            xy = xy(1:4);
            entrada = xy;
            X = [ entrada.* V_max; 1]; % bias
            s1 = W1 * X;
            y1=[pmntanh(s1); 1];
            s2 = W2 * y1;
            y2 = s2;
            Q(acc)=indice_g_motor(x_ini',uf(acc))+y2;
        end
        [val lugar]=min(Q(:));
        Ya(posi_x)=uf(lugar);
        PHI(:,posi_x)=x_ini'.* V_max;
    end
    [W1a,W2a,PI_vector,iteration,lambda]=marq(NetDef,W1a,W2a,PHI,Ya',trparms);
    sal(1)=CI;
    costo=0;
    entrada = [0,0,0,CI]';
    for k=1:etapas-1
        X = [ entrada.* V_max; 1]; % bias
        s1 = W1a * X;
        y1=[pmntanh(s1); 1];
        s2 = W2a * y1;
        xy=mopdm2_motor(t_etapa,entrada,torque,s2,Vf);
        xy = xy(1:4);
        costo=costo+indice_g_motor(entrada,s2);
        entrada = xy;
    end
    evoluc(iterac)=costo; costo
end
Evolc_etapas_motor;
VX;
% figure;
% semilogy(evoluc,'.k');
% xlabel('Iteraciones');
% title('Costo para ir desde x(0)=2 a x(5)=1');
% costo=zeros(4,etapas);estado=costo;u_opt=costo;
