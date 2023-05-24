% Programación dinámica.
% Apredizaje Q. Para péndulo inestable y estable
% con el funcional de costo J=sum((2+u(k))*exp(-x(k)));
%Autor JAP
clear,clc,close all;

% PENDULO EN EL EQUILIBRIO INESTABLE

% Condidiones iniciales
% x = [delta(i); delta_p(i); fi(i); fi_p(i)]
x_ini = [0;0;.2;0];
angulo(1) = x_ini(3);
vfinal=0;

% Algoritmo aprendizaje Q
Max_it=10; % Máximas iteraciones
TM=500; % Cantidad de estados
Nmax=10; % Valores que tiene x
Mmax=15; % Valores que tiene uf
t_etapa = 0.1;
color='.-k';
tic;dx=Nmax; du=Mmax; etapas=10;umin=-100; umax=100;xmin=0; xmax=3;
%%%Carga de datos
rand('state',0);
randn('state',0);

%VALORES OBTENIDOS CON EL PLANO DE FASE
delta= 1*(randn(TM,1));
delta_p = 1*(randn(TM,1));
fi = 2*(randn(TM,1));
fi_p = 2*(randn(TM,1));

M_est = [delta, delta_p, fi, fi_p]; % Matriz de estados

CI=x_ini(3); longx = 4; dd = 1; NP = dx*etapas;
ho = 7;
dW1 = zeros(ho,longx+1); dW2 = zeros(dd,ho+1);
W1 = rand(ho,longx+1); W2 = rand(dd,ho+1); W1a=W1-.5;W2a=W2-.5;
NetDef = ['HHHHHHH'
    'L------'];
maxiter = 10; stop_crit = 1e-4; lambda=1; D=0.01;
trparms=[maxiter stop_crit lambda D]; vfinal=1;
gama(1)=1;
tic; divi=1; xx=0:etapas-1;

Ax=(xmax-xmin)/(Nmax-1); %Delta x
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
        x_ini = M_est(posi_x,:);
        entrada = x_ini;
        X = [ entrada'; 1]; % bias
        s1 = W1a * X;
        y1=[pmntanh(s1); 1];
        s2 = W2a * y1;
        u = s2;
        xy=mopdm2(t_etapa,entrada,u);
        X = [xy; 1];
        s1 = W1 * X;
        y1=[pmntanh(s1); 1];
        s2 = W2 * y1;
        y2=s2;
        C(posi_x)=y2+indice_g(entrada',u);
        entrada = xy;
        PHI(:,posi_x)=x_ini;
        Yo(posi_x)=C(posi_x);
    end
    [W1,W2,PI_vector,iteration,lambda]=marq(NetDef,W1,W2,PHI,Yo,trparms);
    
    for posi_x=1:TM
        x_ini = M_est(posi_x,:);
        entrada = x_ini;
        for acc=1:du
            xy=mopdm2(t_etapa,x_ini,uf(acc));
            entrada = xy;
            X = [ entrada; 1]; % bias
            s1 = W1 * X;
            y1=[pmntanh(s1); 1];
            s2 = W2 * y1;
            y2 = s2;
            Q(acc)=indice_g(x_ini',uf(acc))+y2;
        end
        [val lugar]=min(Q(:));
        J(posi_x)=val;
        Yo(posi_x)=uf(lugar);
        PHI(:,posi_x)=x_ini;
    end
    [W1a,W2a,PI_vector,iteration,lambda]=marq(NetDef,W1a,W2a,PHI,Yo,trparms);
    sal(1)=CI;
    costo=0;
    entrada = [0,0,.1,0]';
    for k=1:etapas-1
        X = [ entrada; 1]; % bias
        s1 = W1a * X;
        y1=[pmntanh(s1); 1];
        s2 = W2a * y1;
        xy=mopdm2(t_etapa,entrada,s2);
        costo=costo+indice_g(entrada,s2);
        entrada = xy;
    end
    evoluc(iterac)=costo; costo
end
Evolucion_etapas;

figure;
semilogy(evoluc,'.k');
xlabel('Iteraciones');
title('Costo para ir desde x(0)=2 a x(5)=1');
costo=zeros(4,etapas);estado=costo;u_opt=costo;
