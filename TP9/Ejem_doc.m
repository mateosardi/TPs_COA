clear,close all,clc;
alfa=3; Nmax=4*alfa; Mmax=3*alfa; color='.-k';
rand('state',0);
tic; dx=Nmax; du=Mmax; etapas=6; xmin=0; xmax=3;umin=-1; umax=1;
% muh=zeros(etapas-1,dx);
Ax=(xmax-xmin)/(Nmax-1); Au=(umax-umin)/(Mmax-1);
for i=1:Nmax
    x(i)=xmin+Ax*(i-1);
end
for i=1:Mmax
    uf(i)=umin+Au*(i-1);
end
CI=2; longx = 2; dd = 1; NP = dx*etapas;
ho = 10;
dW1 = zeros(ho,longx+1); dW2 = zeros(dd,ho+1);
W1 = rand(ho,longx+1); W2 = rand(dd,ho+1); W1a=W1-.5;W2a=W2-.5;
NetDef = ['HHHHHHHHHH'
    'L---------'];
maxiter = 15; stop_crit = 1e-9; lambda=1; D=0;
trparms=[maxiter stop_crit lambda D]; vfinal=1;
gama(1)=1;
tic; divi=1; xx=0:etapas-1;
for iterac=1:100
    Q = 10e3*ones(etapas,dx,du);
    m=zeros(etapas,dx);
    C=zeros(etapas,dx);
    J=zeros(etapas,dx);
    pos=0;
    for posi_x=1:dx
        for k=1:etapas
            if k<etapas
                entrada = [k; x(posi_x)];
                X = [ entrada; 1]; % bias
                s1 = W1a * X;
                y1=[pmntanh(s1); 1];
                s2 = W2a * y1;
                u = s2;
                xy=mopdm(k,x(posi_x),u);
                if (xy>=xmin)&(xy<=xmax)
                    m(k,posi_x)=m(k,posi_x)+1;
                    entrada = [k+1; xy];
                    X = [ entrada; 1];
                    s1 = W1 * X;
                    y1=[pmntanh(s1); 1];
                    s2 = W2 * y1;
                    y2=s2;
                    C(k,posi_x)=y2+indice(k,x(posi_x),u);
                end
            else
                m(k,posi_x)=m(k,posi_x)+1;
                C(k,posi_x)=abs(x(posi_x)-vfinal);
            end
            pos=pos+1;
            PHI(1,pos)=k;
            PHI(2,pos)=x(posi_x);
            Yo(pos)=C(k,posi_x);
        end
    end
    [W1,W2,PI_vector,iteration,lambda]=marq(NetDef,W1,W2,PHI,Yo,trparms);
    for k=1:etapas-1
        for i=1:dx
            for acc=1:du
                xy=mopdm(k,x(i),uf(acc));
                if (xy>=xmin)&(xy<=xmax)
                    [val lugar]=min(abs(x-xy));
                    entrada = [k+1; xy];
                    X = [ entrada; 1]; % bias
                    s1 = W1 * X;
                    y1=[pmntanh(s1); 1];
                    s2 = W2 * y1;
                    y2 = s2;
                    Q(k,i,acc)=indice(k,x(i),uf(acc))+y2;
                end
            end
        end
    end
    pos=0; PHI=0; Yo=0;
    for i=1:dx
        for k=1:etapas-1
            [val lugar]=min(Q(k,i,:));
            J(k,i)=val;
            pos=pos+1;
            Yo(pos)=uf(lugar);
            PHI(1,pos)=k;
            PHI(2,pos)=x(i);
        end
    end
    [W1a,W2a,PI_vector,iteration,lambda]=marq(NetDef,W1a,W2a,PHI,Yo,trparms);
    gama(iterac+1)=100/(100+1.5*iterac);
    sal(1)=CI;
    costo=0;
    for k=1:etapas-1
        entrada = [k; sal(k)];
        X = [ entrada; 1]; % bias
        s1 = W1a * X;
        y1=[pmntanh(s1); 1];
        s2 = W2a * y1;
        consigna(k)=s2;
        sal(k+1)=mopdm(k,sal(k),consigna(k));
        costo(k+1)=costo(k)+indice(k,sal(k),consigna(k));
    end
    costo(k+1)=costo(k+1)+abs(sal(k+1)-vfinal);
    evoluc(iterac)=costo(etapas); costo(etapas)
end
toc
costo=zeros(4,etapas);estado=costo;u_opt=costo;
for CI=1:4
    in=CI-1; estado(CI,1)=in;
    for k=1:etapas-1
        entrada = [k; in];
        X = [ entrada; 1]; % bias
        s1 = W1a * X;
        y1=[pmntanh(s1); 1];
        an = W2a * y1;
        u_opt(CI,k)=an;
        costo(CI,k+1)=indice(k,in,an)+costo(CI,k);
        in=mopdm(k,in,an);
        estado(CI,k+1)=in;
        
        
        VX(:,k)=sx;
       entrada = [sx(1); sx(2);sx(3);sx(4)]+0*diag([sigma_p 0 sigma_tita 0])*randn(4,1);
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
    costo(CI,k+1)=costo(CI,k+1)+abs(in-1);
end
figure;
xx=0:etapas-1;
color='.-';
subplot(2,2,1),plot(xx,estado,color),title('Estados'),ylabel('Estados'),grid on;
subplot(2,2,2),plot(xx,costo,color),title('Costo'),ylabel('Costos'),grid on;
subplot(2,2,3),plot(xx(1:length(u_opt(1,:))),u_opt,color),title('Acci?e control'),ylabel('Acciones de control'),xlabel('Etapas'),grid on,hold on,axis([0 5.1 -1.1 1.1]);
subplot(2,2,4),semilogy(evoluc),title('Evolucion del J(0)');xlabel('Iteraciones');
