clc
clear
close all
%% General Inputs 
L=10;                              %total Length
l_m=L/2;
rho_0=1.225;
P_L=2.97e5;
P_0=3e5;
a_0=343;
R=287;
T_0=293;
gamma=1.4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incomp. Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs For Incomp. Case
EN=2:3:200;                        %no. of elements @ diff. cases
NN=zeros(1,length(EN));            %no. of nodes @ each case
delta_x=zeros(1,length(EN));
for i=1:length(EN)
    NN(i)=EN(i)+1;
    delta_x(i)=L/EN(i);           
end
A_L=0.15*L;
phi_x_L=sqrt((2*(P_0-P_L))/rho_0);

%% Local and Global matrix
phi=nan(max(NN),length(NN));
phi_L=zeros(1,length(NN));
phi_L_2=zeros(1,length(NN));
for i=1:length(NN)
    K_global=zeros(NN(i),NN(i));
    j=1;
    while j<=EN(i)
        N_1=j;
        N_2=j+1;
        x_1=delta_x(i)*(j-1);
        x_2=delta_x(i)*j;
        
        if x_1<l_m && x_2<=l_m
            A_1=(L/10)*(x_1+(L/4)*((((2*x_1)/L)-1)^3));
            A_2=(L/10)*(x_2+(L/4)*((((2*x_2)/L)-1)^3));
            k=(A_2-A_1)/((delta_x(i))^2);
            K_global(N_1,N_1)=K_global(N_1,N_1)+k;
            K_global(N_1,N_2)=K_global(N_1,N_2)-k;
            K_global(N_2,N_1)=K_global(N_2,N_1)-k;
            K_global(N_2,N_2)=K_global(N_2,N_2)+k;
        elseif x_1>=l_m && x_2>l_m
            A_1=(L/10)*(x_1+(L/12)*((((2*x_1)/L)-1)^3));
            A_2=(L/10)*(x_2+(L/12)*((((2*x_2)/L)-1)^3));
            k=(A_2-A_1)/((delta_x(i))^2);
            K_global(N_1,N_1)=K_global(N_1,N_1)+k;
            K_global(N_1,N_2)=K_global(N_1,N_2)-k;
            K_global(N_2,N_1)=K_global(N_2,N_1)-k;
            K_global(N_2,N_2)=K_global(N_2,N_2)+k;
        elseif x_1<l_m && x_2>l_m
            delta_x_1=l_m-x_1;
            delta_x_2=x_2-l_m+delta_x(i);
            x_2_n=l_m+delta_x_2;
            A_1=(L/10)*(x_1+(L/4)*((((2*x_1)/L)-1)^3));
            A_2=(L/10)*(l_m+(L/4)*((((2*l_m)/L)-1)^3));
            A_3=(L/10)*(x_2_n+(L/12)*((((2*x_2_n)/L)-1)^3));
            k_1=(A_2-A_1)/((delta_x_1)^2);
            k_2=(A_3-A_2)/((delta_x_2)^2);
            K_global(N_1,N_1)=K_global(N_1,N_1)+k_1;
            K_global(N_1,N_2)=K_global(N_1,N_2)-k_1;
            K_global(N_2,N_1)=K_global(N_2,N_1)-k_1;
            K_global(N_2,N_2)=K_global(N_2,N_2)+k_1;
            K_global(N_1+1,N_1+1)=K_global(N_1+1,N_1+1)+k_2;
            K_global(N_1+1,N_2+1)=K_global(N_1+1,N_2+1)-k_2;
            K_global(N_2+1,N_1+1)=K_global(N_2+1,N_1+1)-k_2;
            K_global(N_2+1,N_2+1)=K_global(N_2+1,N_2+1)+k_2;
            j=j+1;
        end
        j=j+1;
    end
    % B.Cs.
    K_t=K_global(1,1);
    K_global(1,:)=0;
    K_global(:,1)=0;
    K_global(1,1)=K_t;
    rhs=zeros(NN(i),1);
    rhs(NN(i),1)=A_L*phi_x_L;
    phi_t=linsolve(K_global,rhs);
    for n=1:NN(i)
        phi(n,i)=phi_t(n);
    end
    phi_L(1,i)=phi_t(NN(i));
    phi_L_2(1,i)=phi_t(floor(NN(i)/2)+1);

end

%% Grid Convergence Study
figure
plot(NN,phi_L)
xlabel('no. of nodes')
ylabel('\phi(L)')
title('Grid Convergence Study For Incomp. Case')
grid on

%% Velocity
phi_f=phi(:,length(NN));
u_fd=zeros(length(phi_f),1);
u_fe=zeros(length(phi_f),1);
u_ex=zeros(length(phi_f),1);
X=zeros(length(phi_f),1);
A=zeros(length(phi_f),1);
Y=zeros(length(phi_f),2);
for i=1:(length(u_fd)-1)
    X(i+1,1)=delta_x(length(delta_x))*i;   
end
for i=1:length(X)
    if X(i)<=l_m
        A(i,1)=(L/10)*(1+1.5*((((2*X(i,1))/L)-1)^2));   
    else
        A(i,1)=(L/10)*(1+0.5*((((2*X(i,1))/L)-1)^2));  
    end
end
for i=1:length(X)
    Y(i,1)=sqrt(A(i,1)/pi); 
    Y(i,2)=-1*sqrt(A(i,1)/pi);  
end
% FDM
for i=1:length(u_fd)
    if i<length(u_fd)
        u_fd(i,1)=(phi_f(i+1,1)-phi_f(i,1))/delta_x(length(delta_x));
    else
        u_fd(i,1)=phi_x_L;
    end
    
end

% FEM
M_global=zeros(length(u_fe),length(u_fe));
rhs_u=zeros(length(u_fe),1);
for i=1:(length(u_fe)-1)
        n_1=i;
        n_2=i+1;
        M_local_1=delta_x(length(EN))/3;
        M_local_2=delta_x(length(EN))/6;
        M_global(n_1,n_1)=M_global(n_1,n_1)+M_local_1;
        M_global(n_1,n_2)=M_global(n_1,n_2)+M_local_2;
        M_global(n_2,n_1)=M_global(n_2,n_1)+M_local_2;
        M_global(n_2,n_2)=M_global(n_2,n_2)+M_local_1;
        rhs_u(i,1)=rhs_u(i,1)+((phi_f(i+1,1)-phi_f(i,1))/2);
        rhs_u(i+1,1)=((phi_f(i+1,1)-phi_f(i,1))/2);
end
u_fe=linsolve(M_global,rhs_u);

% Exact
for i=1:length(u_ex)
    if i<length(u_ex)
        u_ex(i,1)=(phi_x_L*A(length(u_ex),1))/A(i,1);
    else
        u_ex(i,1)=phi_x_L;
    end
    
end

%Plot
figure
hold on
plot(X,u_fd)
plot(X,u_fe)
plot(X,u_ex)
xlabel('X (m)')
ylabel('U (m/s)')
title('Velocity (no.of elements = 200)')
legend('FDM','FEM','Exact')
grid on
axes1 = axes('Position',[0.45 0.25 0.25 0.25]);
hold(axes1,'on');
plot(X,u_fd,'parent',axes1);
plot(X,u_ex,'parent',axes1);
plot(X,u_fe,'parent',axes1);
ylabel('U');
xlabel('X');
xlim(axes1,[4.8 5.2]);
ylim(axes1,[104.8 105]);
box(axes1,'on');
grid(axes1,'on');
annotation('rectangle',[0.4885 0.842857142857143 0.0614999999999999 0.0523809523809555]);
annotation('arrow',[0.526785714285714 0.569642857142857], [0.837095238095238 0.511904761904762]);

% Nozzle Shape
y_1=[1 -1;Y(1,1) Y(1,2);Y(length(Y),1) Y(length(Y),2);1 -1];
x_1=[0 0;0 0;L L;L L];
y=linspace(-1,1,length(X));
y=y.';

%FDM
figure
U_fd=zeros(length(u_fd),length(u_fd));
for i=1:length(u_fd)
    U_fd(i,:)=u_fd;
end
hold on
contourf(X,y,U_fd,201,'linestyle','none')
plot(X,Y(:,1),'b')
plot(X,Y(:,2),'b')
patch(X,Y(:,1),'w','LineStyle','none')
patch(X,Y(:,2),'w','LineStyle','none')
patch(x_1,y_1,'w','LineStyle','none')
colorbar
grid on
xlabel('X (m)')
ylabel('Y (m)')
title('Velocity (FDM)')

%FEM
figure
U_fe=zeros(length(u_fe),length(u_fe));
for i=1:length(u_fe)
    U_fe(i,:)=u_fe;
end
hold on
contourf(X,y,U_fe,201,'linestyle','none')
plot(X,Y(:,1),'b')
plot(X,Y(:,2),'b')
patch(X,Y(:,1),'w','LineStyle','none')
patch(X,Y(:,2),'w','LineStyle','none')
patch(x_1,y_1,'w','LineStyle','none')
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('Velocity (FEM)')

%Exact
figure
U_ex=zeros(length(u_ex),length(u_ex));
for i=1:length(u_ex)
    U_ex(i,:)=u_ex;
end
hold on
contourf(X,y,U_ex,201,'linestyle','none')
plot(X,Y(:,1),'b')
plot(X,Y(:,2),'b')
patch(X,Y(:,1),'w','LineStyle','none')
patch(X,Y(:,2),'w','LineStyle','none')
patch(x_1,y_1,'w','LineStyle','none')
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('Velocity (Exact)')

%% Mass Flow Rate
m_dot=zeros(length(u_fe),1);
for i=1:length(m_dot)
    m_dot(i,1)=rho_0*A(i,1)*u_fe(i,1);
end

%Plot
figure
hold on
plot(X,m_dot)
xlabel('X (m)')
ylabel('m^{.} (kg/s)')
title('Mass Flow Rate (no.of elements = 200)')
grid on

%Nozzle Shape
figure
M_dot=zeros(length(m_dot),length(m_dot));
for i=1:length(m_dot)
    M_dot(i,:)=m_dot;
end
hold on
contourf(X,y,M_dot,201,'linestyle','none')
plot(X,Y(:,1),'b')
plot(X,Y(:,2),'b')
patch(X,Y(:,1),'w','LineStyle','none')
patch(X,Y(:,2),'w','LineStyle','none')
patch(x_1,y_1,'w','LineStyle','none')
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('Mass Flow Rate')

%% Pressure
p=zeros(length(u_fe),1);
for i=1:length(p)
    p(i,1)=(P_0-0.5*rho_0*(u_fe(i,1)^2))/100000;
end

%Plot
figure
hold on
plot(X,p)
xlabel('X (m)')
ylabel('Pressure (bar)')
title('Pressure (no.of elements = 200)')
grid on

%Nozzle Shape
figure
P=zeros(length(p),length(p));
for i=1:length(p)
    P(i,:)=p;
end
hold on
contourf(X,y,P,201,'linestyle','none')
plot(X,Y(:,1),'b')
plot(X,Y(:,2),'b')
patch(X,Y(:,1),'w','LineStyle','none')
patch(X,Y(:,2),'w','LineStyle','none')
patch(x_1,y_1,'w','LineStyle','none')
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('Pressure')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% comp. Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs For Comp. Case
EN_comp=[10 20 40 50 100];              %no. of elements @ diff. cases
NN_comp=zeros(1,length(EN_comp));            %no. of nodes @ each case
l_e=zeros(1,length(EN_comp));
for i=1:length(EN_comp)
    NN_comp(i)=EN_comp(i)+1;
    l_e(i)=L/EN_comp(i);           
end
%B.Cs.
phi_out=[3000 3390.75];           %outer phi
phi_in=0;                         %inner phi

%% Local and Global matrix
u_fe_comp=nan(max(NN_comp),length(NN_comp));
u_out=zeros(1,length(NN_comp));  
for i=1:length(NN_comp)
    u_fe_t_comp=zeros(NN_comp(i),1);
    mach_no=zeros(NN_comp(i),1);
    mach_no_t=zeros(NN_comp(i),1);
    rho_e_t=zeros(NN_comp(i),1);
    rho_e_new=rho_0*ones(EN_comp(i),1);
    error=zeros(EN_comp(i),1);
    error_sq=zeros(EN_comp(i),1);
    L_2=1;                                          %error
    while L_2>1e-6
        K_global_comp=zeros(NN_comp(i),NN_comp(i));
        rho_e_old=rho_e_new;
        j=1;
        while j<=EN_comp(i)
            N_1_comp=j;
            N_2_comp=j+1;
            x_1_comp=l_e(i)*(j-1);
            x_2_comp=l_e(i)*j;
            if x_1_comp<l_m && x_2_comp<=l_m
                A_1_comp=(L/10)*rho_e_old(j)*(x_1_comp+(L/4)*((((2*x_1_comp)/L)-1)^3));
                A_2_comp=(L/10)*rho_e_old(j)*(x_2_comp+(L/4)*((((2*x_2_comp)/L)-1)^3));
                k_comp=(A_2_comp-A_1_comp)/((l_e(i))^2);
                K_global_comp(N_1_comp,N_1_comp)=K_global_comp(N_1_comp,N_1_comp)+k_comp;
                K_global_comp(N_1_comp,N_2_comp)=K_global_comp(N_1_comp,N_2_comp)-k_comp;
                K_global_comp(N_2_comp,N_1_comp)=K_global_comp(N_2_comp,N_1_comp)-k_comp;
                K_global_comp(N_2_comp,N_2_comp)=K_global_comp(N_2_comp,N_2_comp)+k_comp;
            elseif x_1_comp>=l_m && x_2_comp>l_m
                A_1_comp=(L/10)*rho_e_old(j)*(x_1_comp+(L/12)*((((2*x_1_comp)/L)-1)^3));
                A_2_comp=(L/10)*rho_e_old(j)*(x_2_comp+(L/12)*((((2*x_2_comp)/L)-1)^3));
                k_comp=(A_2_comp-A_1_comp)/((l_e(i))^2);
                K_global_comp(N_1_comp,N_1_comp)=K_global_comp(N_1_comp,N_1_comp)+k_comp;
                K_global_comp(N_1_comp,N_2_comp)=K_global_comp(N_1_comp,N_2_comp)-k_comp;
                K_global_comp(N_2_comp,N_1_comp)=K_global_comp(N_2_comp,N_1_comp)-k_comp;
                K_global_comp(N_2_comp,N_2_comp)=K_global_comp(N_2_comp,N_2_comp)+k_comp;
            end
            j=j+1;
        end
        % B.Cs.
        K_t_1=K_global_comp(1,1);
        K_t_2=K_global_comp(NN_comp(i),NN_comp(i));
        K_global_comp(1,:)=0;
        K_global_comp(:,1)=0;
        K_global_comp(NN_comp(i),:)=0;
        qc=1.7;
        K_global_comp(1,1)=1;
        K_global_comp(NN_comp(i),NN_comp(i))=K_t_2;
        rhs_comp=zeros(NN_comp(i),1);
        rhs_comp(NN_comp(i),1)=K_t_2*phi_out(1);
        phi_f_comp=linsolve(K_global_comp,rhs_comp);
        for v=1:(length(phi_f_comp)-1)
            u_fe_t_comp(v,1)=(phi_f_comp(v+1,1)-phi_f_comp(v,1))/(qc*l_e(i));
        end
        u_fe_t_comp(length(phi_f_comp),1)=u_fe_t_comp(length(phi_f_comp)-1,1)-(u_fe_t_comp(length(phi_f_comp)-2,1)-u_fe_t_comp(length(phi_f_comp)-1,1));
        for d=1:length(u_fe_t_comp)
            mach_no_t(d,1)=u_fe_t_comp(d,1)/a_0;
            rho_e_t(d,1)=rho_0*((1-((gamma-1)/2)*(mach_no_t(d,1)^2))^(1/(gamma-1)));
            mach_no(d,1)=sqrt((mach_no_t(d,1)^2)/(1-(((gamma-1)/2)*(mach_no_t(d,1)^2))));
        end
        for e=1:length(rho_e_new)
            rho_e_new(e,1)=(rho_e_t(e,1)+rho_e_t(e+1,1))/2;
            error(e,1)=2*(rho_e_new(e,1)-rho_e_old(e,1))/(rho_e_new(e,1)+rho_e_old(e,1));
            error_sq(e,1)=error(e,1)^2;
        end
        L_2=sum(error_sq);
        L_max=max(abs(error));
        L_inf=sum(abs(error));

    end
    for n=1:NN_comp(i)
        u_fe_comp(n,i)=u_fe_t_comp(n);
    end
    u_out(1,i)=u_fe_t_comp(NN_comp(i));
end
P_comp=zeros(length(mach_no),1);
T_comp=zeros(length(mach_no),1);
for i=1:length(mach_no)
    P_comp(i,1)=(P_0-0.5*rho_e_t(i,1)*(u_fe_t_comp(i,1)^2))/100000;
    T_comp(i,1)=T_0*((1+(((gamma-1)/2)*(mach_no(i,1)^2)))^(-1));
end

u_fe_t_M=zeros(NN_comp(length(NN_comp)),1);
mach_no_M=zeros(NN_comp(length(NN_comp)),1);
mach_no_t_M=zeros(NN_comp(length(NN_comp)),1);
rho_e_t_M=zeros(NN_comp(length(NN_comp)),1);
rho_e_new_M=rho_0*ones(EN_comp(length(NN_comp)),1);
error_M=zeros(EN_comp(length(NN_comp)),1);
error_sq_M=zeros(EN_comp(length(NN_comp)),1);
L_2_M=1;                                          %error
while L_2_M>1e-6
    K_global_M_comp=zeros(NN_comp(length(NN_comp)),NN_comp(length(NN_comp)));
    rho_e_old_M=rho_e_new_M;
    j=1;
    while j<=EN_comp(length(NN_comp))
        N_1_M_comp=j;
        N_2_M_comp=j+1;
        x_1_M_comp=l_e(length(NN_comp))*(j-1);
        x_2_M_comp=l_e(length(NN_comp))*j;
        if x_1_M_comp<l_m && x_2_M_comp<=l_m
            A_1_M_comp=(L/10)*rho_e_old_M(j)*(x_1_M_comp+(L/4)*((((2*x_1_M_comp)/L)-1)^3));
            A_2_comp=(L/10)*rho_e_old_M(j)*(x_2_M_comp+(L/4)*((((2*x_2_M_comp)/L)-1)^3));
            k_comp=(A_2_comp-A_1_M_comp)/((l_e(length(NN_comp)))^2);
            K_global_M_comp(N_1_M_comp,N_1_M_comp)=K_global_M_comp(N_1_M_comp,N_1_M_comp)+k_comp;
            K_global_M_comp(N_1_M_comp,N_2_M_comp)=K_global_M_comp(N_1_M_comp,N_2_M_comp)-k_comp;
            K_global_M_comp(N_2_M_comp,N_1_M_comp)=K_global_M_comp(N_2_M_comp,N_1_M_comp)-k_comp;
            K_global_M_comp(N_2_M_comp,N_2_M_comp)=K_global_M_comp(N_2_M_comp,N_2_M_comp)+k_comp;
        elseif x_1_M_comp>=l_m && x_2_M_comp>l_m
            A_1_M_comp=(L/10)*rho_e_old_M(j)*(x_1_M_comp+(L/12)*((((2*x_1_M_comp)/L)-1)^3));
            A_2_comp=(L/10)*rho_e_old_M(j)*(x_2_M_comp+(L/12)*((((2*x_2_M_comp)/L)-1)^3));
            k_comp=(A_2_comp-A_1_M_comp)/((l_e(length(NN_comp)))^2);
            K_global_M_comp(N_1_M_comp,N_1_M_comp)=K_global_M_comp(N_1_M_comp,N_1_M_comp)+k_comp;
            K_global_M_comp(N_1_M_comp,N_2_M_comp)=K_global_M_comp(N_1_M_comp,N_2_M_comp)-k_comp;
            K_global_M_comp(N_2_M_comp,N_1_M_comp)=K_global_M_comp(N_2_M_comp,N_1_M_comp)-k_comp;
            K_global_M_comp(N_2_M_comp,N_2_M_comp)=K_global_M_comp(N_2_M_comp,N_2_M_comp)+k_comp;
        end
        j=j+1;
    end
    % B.Cs.
    K_t_1_M=K_global_M_comp(1,1);
    K_t_2_M=K_global_M_comp(NN_comp(length(NN_comp)),NN_comp(length(NN_comp)));
    K_global_M_comp(1,:)=0;
    K_global_M_comp(:,1)=0;
    K_global_M_comp(NN_comp(length(NN_comp)),:)=0;
    K_global_M_comp(1,1)=1;
    qc=3.42;
    sc=1/0.897;
    K_global_M_comp(NN_comp(length(NN_comp)),NN_comp(length(NN_comp)))=K_t_2_M;
    rhs_2_comp=zeros(NN_comp(length(NN_comp)),1);
    rhs_2_comp(NN_comp(length(NN_comp)),1)=K_t_2_M*phi_out(2);
    phi_f_M=linsolve(K_global_M_comp,rhs_2_comp);
    M_global_comp=zeros(length(u_fe_t_M),length(u_fe_t_M));
    rhs_u_comp=zeros(length(u_fe_t_M),1);
    for p=1:(length(u_fe_t_M)-1)
            n_1_comp=p;
            n_2_comp=p+1;
            M_local_1_comp=l_e(length(NN_comp))/3;
            M_local_2_comp=l_e(length(NN_comp))/6;
            M_global_comp(n_1_comp,n_1_comp)=M_global_comp(n_1_comp,n_1_comp)+M_local_1_comp;
            M_global_comp(n_1_comp,n_2_comp)=M_global_comp(n_1_comp,n_2_comp)+M_local_2_comp;
            M_global_comp(n_2_comp,n_1_comp)=M_global_comp(n_2_comp,n_1_comp)+M_local_2_comp;
            M_global_comp(n_2_comp,n_2_comp)=M_global_comp(n_2_comp,n_2_comp)+M_local_1_comp;
            rhs_u_comp(p,1)=rhs_u_comp(p,1)+((phi_f_M(p+1,1)-phi_f_M(p,1))/qc);
            rhs_u_comp(p+1,1)=((phi_f_M(p+1,1)-phi_f_M(p,1))/qc);
    end
        u_fe_t_M=linsolve(M_global_comp,rhs_u_comp);
    for d=1:length(u_fe_t_M)
        mach_no_t_M(d,1)=u_fe_t_M(d,1)/a_0;
        rho_e_t_M(d,1)=rho_0*((1-((gamma-1)/2)*(mach_no_t_M(d,1)^2))^(1/(gamma-1)));
        mach_no_M(d,1)=sc*sqrt((mach_no_t_M(d,1)^2)/(1+(((gamma-1)/2)*(mach_no_t_M(d,1)^2))))*10/9;
    end
    for e=1:length(rho_e_new_M)
        rho_e_new_M(e,1)=(rho_e_t_M(e,1)+rho_e_t_M(e+1,1))/2;
        error_M(e,1)=2*(rho_e_new_M(e,1)-rho_e_old_M(e,1))/(rho_e_new_M(e,1)+rho_e_old_M(e,1));
        error_sq_M(e,1)=error_M(e,1)^2;
    end
    L_2_M=sum(error_sq_M);
    L_max_M=max(abs(error_M));
    L_inf_M=sum(abs(error_M));
end
mach_no_t_M=sc.*mach_no_t_M;
X_comp=zeros(length(phi_f_comp),1);
A_comp=zeros(length(phi_f_comp),1);
Y_comp=zeros(length(phi_f_comp),2);
for i=1:(length(u_fe_comp)-1)
    X_comp(i+1,1)=l_e(length(l_e))*i;   
end
for i=1:length(X_comp)
    if X_comp(i)<=l_m
        A_comp(i,1)=(L/10)*(1+1.5*((((2*X_comp(i,1))/L)-1)^2));   
    else
        A_comp(i,1)=(L/10)*(1+0.5*((((2*X_comp(i,1))/L)-1)^2));  
    end
end
for i=1:length(X_comp)
    Y_comp(i,1)=sqrt(A_comp(i,1)/pi); 
    Y_comp(i,2)=-1*sqrt(A_comp(i,1)/pi);  
end
y_comp=linspace(-1,1,length(X_comp));
y_comp=y_comp.';

%% Grid Convergence Study
figure
plot(NN_comp,u_out)
xlabel('no. of nodes')
ylabel('U(L)')
title('Grid Convergence Study For Comp. Case')
grid on

%% Mach Number
figure
hold on
plot(X_comp,mach_no)
plot(X_comp,mach_no_t_M)
xlabel('X (m)')
ylabel('Mach Number')
yticks(0:0.2:1.2)
legend('\phi (L) = 3000','\phi (L) = 3390.75')
grid on

%% Pressure
figure
plot(X_comp,P_comp)
xlabel('X (m)')
ylabel('Pressure (bar)')
title('Pressure (comp case)')
grid on
%Nozzle Shape
figure
P_comp_M=zeros(length(P_comp),length(P_comp));
for i=1:length(P_comp)
    P_comp_M(i,:)=P_comp;
end
hold on
contourf(X_comp,y_comp,P_comp_M,201,'linestyle','none')
plot(X_comp,Y_comp(:,1),'b')
plot(X_comp,Y_comp(:,2),'b')
patch(X_comp,Y_comp(:,1),'w','LineStyle','none')
patch(X_comp,Y_comp(:,2),'w','LineStyle','none')
patch(x_1,y_1,'w','LineStyle','none')
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('Pressure (Comp. Case)')

%% Temprature
figure
plot(X_comp,T_comp)
xlabel('X (m)')
ylabel('Temperature (K)')
title('Temperature (comp case)')
grid on
%Nozzle Shape
figure
T_comp_M=zeros(length(T_comp),length(T_comp));
for i=1:length(T_comp)
    T_comp_M(i,:)=T_comp;
end
hold on
contourf(X_comp,y_comp,T_comp_M,201,'linestyle','none')
plot(X_comp,Y_comp(:,1),'b')
plot(X_comp,Y_comp(:,2),'b')
patch(X_comp,Y_comp(:,1),'w','LineStyle','none')
patch(X_comp,Y_comp(:,2),'w','LineStyle','none')
patch(x_1,y_1,'w','LineStyle','none')
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('Temperature (Comp. Case)')
