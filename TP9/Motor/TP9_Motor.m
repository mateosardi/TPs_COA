% Programación dinámica.
%Autor JAP
clear,clc,close all;
% MOTOR DE CC
Pr=[   0.4127    0.0103    0.1830
    0.0103    0.0012    0.0219
    0.1830    0.0219   17.3021]/1000;


% Condidiones iniciales
% x = [ia(i); w(i); theta(i)]
x_ini = [0;0;2];
angulo(1) = x_ini(3);

% Algoritmo aprendizaje Q
Max_it=6; % Máximas iteraciones
TM=100; % Cantidad de estados
Mmax=3; % Valores que tiene uf
t_etapa = 1e-4;
color='.-k';
tic; du=Mmax;
etapas=50;umin=-10; umax=10;
%%%Carga de datos
rand('state',0);
randn('state',0);

%%%Carga de datos
rand('state',0);
randn('state',0);

%VALORES OBTENIDOS CON EL PLANO DE FASE
% Valores máximos
ia_m = 20e-3;
w_m = 40;
theta_m = 2;
V_max = [1/ia_m; 1/w_m; 1/theta_m]; % Vector de valores máximos
% V_max = [1; 1; 1];

ia= ia_m*(randn(TM,1));
w = w_m*(randn(TM,1));
theta = theta_m*(randn(TM,1));

M_est = [ia, w, theta]; % Matriz de estados

CI=x_ini(3); longx = 3; dd = 1;
ho = 7;
dW1 = zeros(ho,longx+1); dW2 = zeros(dd,ho+1);
W1 = rand(ho,longx+1); W2 = rand(dd,ho+1); W1a=W1-.5;W2a=W2-.5;
NetDef = ['HHHHHHH'
    'L------'];
maxiter = 10; stop_crit = 1e-2; lambda=1; D=0;
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
            xy=mopdm2_motor(t_etapa,entrada,u);
            %         X = [xy; 1];
            %         s1 = W1 * X;
            %         y1=[pmntanh(s1); 1];
            %         s2 = W2 * y1;
            %         y2=s2;
            %         C(posi_x)=y2+indice_g_motor(entrada',u);
            C_C=C_C+indice_g_motor(entrada,u);
            entrada = xy;
        end
        %         PHI(:,posi_x)=x_ini;
        Yo(posi_x)=C_C;
    end
    [W1,W2,PI_vector,iteration,lambda]=marq(NetDef,W1,W2,PHI,Yo,trparms);
    C_Yo = Yo;
    Ya = zeros(TM,1);
    for posi_x=1:TM
        x_ini = M_est(posi_x,:);
        entrada = x_ini;
        for acc=1:du
            xy=mopdm2_motor(t_etapa,x_ini,uf(acc));
            entrada = xy;
            X = [ entrada.* V_max; 1]; % bias
            s1 = W1 * X;
            y1=[pmntanh(s1); 1];
            s2 = W2 * y1;
            y2 = s2;
            Q(acc)=indice_g_motor(x_ini',uf(acc))+y2;
            %             Q(acc)=indice_g_motor(x_ini',uf(acc))+xy'*Pr*xy;
        end
        [val lugar]=min(Q(:));
        J(posi_x)=val;
        Ya(posi_x)=uf(lugar);
        PHI(:,posi_x)=x_ini'.* V_max;
    end
    [W1a,W2a,PI_vector,iteration,lambda]=marq(NetDef,W1a,W2a,PHI,Ya',trparms);
    sal(1)=CI;
    costo=0;
    entrada = [0,0,CI]';
    for k=1:etapas-1
        X = [ entrada.* V_max; 1]; % bias
        s1 = W1a * X;
        y1=[pmntanh(s1); 1];
        s2 = W2a * y1;
        xy=mopdm2_motor(t_etapa,entrada,s2);
        costo=costo+indice_g_motor(entrada,s2);
        entrada = xy;
    end
    evoluc(iterac)=costo; costo
end
Evoluc_etapas_motor;

figure;
semilogy(evoluc,'.k');
xlabel('Iteraciones');
title('Costo para ir desde x(0)=2 a x(5)=1');
costo=zeros(4,etapas);estado=costo;u_opt=costo;
