%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Linear SM-MPSC for N-body MSD
%
% Alexandre Didier
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all
rng(0)

show_figures=true;
feasibility_grid_size=0.1;
use_lorenzen=false;
k_only=true;
n_MC=1e5;
use_PGSD=true;
compare_linearMPSC=true;
%% Setup
% Number of MSD
n_MSD=3;
if n_MSD>3
    show_figures=false;
    compare_linearMPSC=false;
end
% Number of states
n=2*n_MSD;

% Number of inputs
m=n_MSD;



% Sample time
T_s=0.2;

% Simulation time
time=30;

% MPSC horizon
N=10;

% Number of parameters
if k_only
    p=n_MSD-1;
else
    p=2*(n_MSD-1);
end

% Parameter bounds
k_min=0.05;
k_max=0.25;
d_min=0.05;
d_max=0.25;

% True parameters
[k_true,d_true]=get_MSD_params(n_MSD, k_min,k_max,d_min,d_max);
if k_only
    d_true=0.1*ones(n_MSD-1,1);
end

% True Dynamics matrices
[A,B]=get_true_dynamics(n_MSD, T_s, k_true, d_true);

% Model Dynamics matrices
A_theta=get_uncertain_dynamics(n_MSD, T_s, k_only, d_true);
B_theta=B;


%Construct Reference
x_ref_PGSD=[];
for i=1:n_MSD
    x_ref_PGSD=[x_ref_PGSD;[2.3*sin(0.1*(1:time/T_s+1)+2*pi*i/n_MSD);zeros(1,time/T_s+1)]];
end

x_init=zeros(n,1);

%% Constraints
% State constraints H_x*x<=h_x
H_x=[];
h_x=[];
for i=1:n_MSD
    H_x=[H_x; zeros(4,2*(i-1)), [1,0;-1,0;0,1;0,-1], zeros(4,2*(n_MSD-i))];
    h_x=[h_x;2.3*ones(4,1)];
end

X=Polyhedron(H_x,h_x);

% Input constraints H_u*u<=h_u
H_u=[];
h_u=[];
for i=1:n_MSD
    H_u=[H_u; zeros(2,i-1), [1;-1], zeros(2,n_MSD-i)];
    h_u=[h_u;3.5*ones(2,1)];
end

U=Polyhedron(H_u,h_u);

%% Set of parameters
% Set of Thetas
H_theta=[];
h_theta=[];
if ~k_only
    for i=1:n_MSD-1
        H_theta=[H_theta; zeros(2,i-1), [1;-1], zeros(2,n_MSD-1-i), zeros(2,n_MSD-1)];
        h_theta=[h_theta;k_max;-k_min];
    end
    for i=1:n_MSD-1
        H_theta=[H_theta; zeros(2,n_MSD-1), zeros(2,i-1), [1;-1], zeros(2,n_MSD-1-i)];
        h_theta=[h_theta;d_max;-d_min];
    end
else
    for i=1:n_MSD-1
        H_theta=[H_theta; zeros(2,i-1), [1;-1], zeros(2,n_MSD-1-i)];
        h_theta=[h_theta;k_max;-k_min];
    end
end

Omega=Polyhedron(H_theta,h_theta);

%% Disturbance set
% Disturbance set

H_w=H_x;
h_w=[0.001*ones(size(H_w,1),1)];
W=Polyhedron(H_w,h_w);



%% Cost matrices
% State cost matrix
Q=diag(repmat([100,30],1,n_MSD))/100;

% Input cost matrix
R=diag(repmat(100,1,n_MSD))/100;

%% Precomputations for SM-MPC according to Köhler et al.
% Centre of set of thetas
for i=1:p
    theta_bar_0(i,1)=(max(Omega.V(:,i))+min(Omega.V(:,i)))/2;
end
theta_hat_0=theta_bar_0;

% Hypercube length
eta_0 = max(Omega.V(:,1))-min(Omega.V(:,1));

% Hypercube vertices
H_hc=H_theta;
h_hc=0.5*ones(size(H_theta,1),1);
HC=Polyhedron(H_hc,h_hc);
e_tilde=HC.V;

% Vertices
vertices_theta=Omega.V;

%% Computation of an invariant Ellipsoid
% Find Feedback K and terminal cost P
X_PK=sdpvar(n,n);
Y_PK=sdpvar(m,n);
lambda_PK=0.995;
objective_PK=-logdet(X_PK);
constraints_PK=[];

for i=1:size(vertices_theta,1)
    A_vertex=A_theta{1};
    if k_only
        for j=1:n_MSD-1
            A_vertex=A_vertex+A_theta{j+1}*vertices_theta(i,j);
        end
    else
        for j=1:n_MSD-1
            A_vertex=A_vertex+A_theta{j+1}*vertices_theta(i,j)+A_theta{j+n_MSD}*vertices_theta(i,j+n_MSD-1);
        end
    end
    % Ricciati equation
    constraints_PK=[constraints_PK,...
        [X_PK,(A_vertex*X_PK+B_theta*Y_PK)', Q^0.5*X_PK, Y_PK'*R^0.5;...
        (A_vertex*X_PK+B_theta*Y_PK),X_PK,zeros(n), zeros(n,m);...
        (Q^0.5*X_PK)',zeros(n),eye(n),zeros(n,m);...
        (Y_PK'*R^0.5)',zeros(m,n),zeros(m,n),eye(m)]>=0];
    %lambda-contractivity
    constraints_PK=[constraints_PK,...
        [lambda_PK*X_PK,(A_vertex*X_PK+B_theta*Y_PK)';...
        (A_vertex*X_PK+B_theta*Y_PK),lambda_PK*X_PK]>=0];
end

% Constraint satisfaction
F=[X.A./X.b;zeros(size(U.A,1),n)];
G=[zeros(size(X.A,1),m);U.A./U.b];
for i=1:size(X.A,1)+size(U.A,1)
    constraints_PK=[constraints_PK,[1,F(i,:)*X_PK+G(i,:)*Y_PK;...
        (F(i,:)*X_PK+G(i,:)*Y_PK)',X_PK]>=0];
end

sol=optimize(constraints_PK,objective_PK,sdpsettings('solver','mosek','verbose',0));
P=inv(value(X_PK));
K=value(Y_PK)*P;

%%  Compute Polytope
% Closed-loop matrices
for i=1:size(vertices_theta,1)
    A_vertex=A_theta{1};
    if k_only
        for j=1:n_MSD-1
            A_vertex=A_vertex+A_theta{j+1}*vertices_theta(i,j);
        end
    else
        for j=1:n_MSD-1
            A_vertex=A_vertex+A_theta{j+1}*vertices_theta(i,j)+A_theta{j+n_MSD}*vertices_theta(i,j+n_MSD-1);
        end
    end
    A_cl{i}=A_vertex+B_theta*K;
end
% Initial set of admissible states
A_s=[X.A./X.b;(U.A./U.b)*K];
b_s=[X.b./X.b;U.b./U.b];

% Some setup
x_0=zeros(2,1);
options=optimoptions(@linprog,'Display','off');
lambda=lambda_PK;
i=1;

% Algorithm to get a lambda-contractive set
while (i<=size(A_s,1))
    %print(i)
    %print(size(A_s,1))
    a=A_s(i,:);
    b=b_s(i,:);
    for j=1:size(vertices_theta,1)
        [~,fval]=linprog((-a*A_cl{j})',A_s,b_s,[],[],[],[],options);
        c_i(j)=-fval-lambda*b;
    end
    for j=1:size(vertices_theta,1)
        if c_i(j)>0
            A_s=[A_s;1/lambda*a*A_cl{j}];
            b_s=[b_s;b];
        end
    end
    i=i+1;
end

% Compute Polyhedron F_0 and define H_f
F_0=Polyhedron(A_s,b_s);
F_0.minHRep();
H_f=F_0.A;
h_f=F_0.b;
if use_lorenzen
    vertices_F=F_0.V;
end

%% More precomputations
% Computation of rho_theta_bar_0
% Note that due to the way F_0 is computed, lambda is an upper bound on
% rho_theta from the paper and the term L_b*eta_0 can be removed.

rho_theta_bar_0=lambda;

% Computation of w_bar
w_bar_i=zeros(size(H_f,1),1);
w0=zeros(n,1);
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon=[];
options=optimoptions(@fmincon,'Display','off');
for i=1:size(H_f,1)
    fun=@(w)-H_f(i,:)*w;
    [~,fval]=fmincon(fun,w0,H_w,h_w,Aeq,beq,lb,ub,nonlcon,options);
    w_bar_i(i)=-fval;
end
w_bar=max(w_bar_i);

% Computation of c_j
c=zeros(size(F,1),1);
x0=zeros(n,1);
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon=[];
options=optimoptions(@fmincon,'Display','off');
for i=1:size(c,1)
    fun=@(x)-(F(i,:)+G(i,:)*K)*x;
    [~,fval]=fmincon(fun,x0,H_f,h_f,Aeq,beq,lb,ub,nonlcon,options);
    c(i)=-fval;
end

% Terminal Set formula if rho_theta_bar_0+max(c)*w_bar<1 (true for this
% example)
if rho_theta_bar_0+max(c)*w_bar>=1
    error('Terminal condition not fulfilled!')
end
H_term=[max(c)*H_f,max(c)*ones(size(H_f,1),1)];
h_term=ones(size(H_term,1),1);



%% Set up optimisation problem using Yalmip
% We use the formulation from
% Lorenzen et al - Robust MPC with recursive model update

% Optimisation Variables
% States and inputs with respect to theta_bar and theta_hat
% Note that z_X are the translations of the polytopes F_i and s the
% dilations
if use_lorenzen
    x_hat=sdpvar(n, 1);
    z_X=sdpvar(n, N+1);
    
    % Input u=Kx+v
    
    u_hat=sdpvar(m, 1);
    v=sdpvar(m, N);
    u_L_k=sdpvar(m,1);
    
    % Dilations
    
    s=sdpvar(1,N+1,'full');
    
    % Parameter
    H_theta_k=sdpvar(2*p,p);
    h_theta_k=sdpvar(2*p,1);
    
    for k=1:N
        for j=1:size(vertices_F,1)
            Delta{k}{j}=sdpvar(size(H_f,1),size(vertices_theta,1));
        end
    end
    
    % Objective function
    
    objective_MPSC=norm(u_L_k(:,1)-u_hat(:,1),2)^2;
    
    
    % Constraints
    % Initial Constraints
    constraints_MPSC=[-H_f*z_X(:,1)-s(:,1)*ones(size(H_f,1),1)<=-H_f*x_hat(:,1)];
    constraints_MPSC=[constraints_MPSC,u_hat(:,1)==v(:,1)+K*x_hat(:,1)];
    % Constraints for every predicted time step
    for k=1:N
        k
        % State and input constraints
        % F and G are state and input constraint matrices so that (x,u)\in
        % {x,u| Fx+Gu<1}
        constraints_MPSC=[constraints_MPSC,(F+G*K)*z_X(:,k)+G*v(:,k)+c*s(:,k)<=ones(size(F,1),1)];
        
        % Tube inclusion constraint
        for j=1:size(vertices_F,1)
            constraints_MPSC=[constraints_MPSC,Delta{k}{j}>=0];
            constraints_MPSC=[constraints_MPSC,Delta{k}{j}*h_theta_k+H_f*(A_theta{1}*(z_X(:,k)+s(:,k)*vertices_F(j,:)')+B_theta*(K*(z_X(:,k)+s(:,k)*vertices_F(j,:)')+v(:,k))-z_X(:,k+1))-s(:,k+1)*ones(size(H_f,1),1)<=-w_bar_i];
            temp_matrix=[];
            if k_only
                for i=1:p
                    temp_matrix=[temp_matrix, A_theta{i+1}*(z_X(:,k)+s(:,k)*vertices_F(j,:)')];
                end
            else
                for i=1:2*(n_MSD-1)
                    temp_matrix=[temp_matrix, A_theta{i+1}*(z_X(:,k)+s(:,k)*vertices_F(j,:)')];
                end
            end
            constraints_MPSC=[constraints_MPSC,H_f*temp_matrix==Delta{k}{j}*H_theta_k];
        end
    end
    
    % Terminal constraint 
    constraints_MPSC=[constraints_MPSC,H_term*[z_X(:,N+1);s(:,N+1)]<=h_term];
    
    % Optimizer
    opt_MPSC=optimizer(constraints_MPSC, objective_MPSC, sdpsettings('solver','mosek'), {u_L_k,x_hat(:,1),H_theta_k, h_theta_k},{u_hat(:,1),z_X,s,x_hat});
else
    objective_MPSC=0;
    % Decision variables
    x_bar=sdpvar(n,N+1,'full');
    v=sdpvar(m,N,'full');
    u_bar=sdpvar(m,N,'full');
    s=sdpvar(1,N+1,'full');
    theta_bar=sdpvar(p,1,'full');
    eta=sdpvar(1,1,'full');
    
    u_L_k=sdpvar(m,1);
    
    
    
    % Quadratic objective
    objective_MPSC=norm(u_L_k(:,1)-u_bar(:,1),2)^2;
    
    % Initial constraint
    constraints_MPSC=[s(1,1)==0];
    
    A_theta_bar=A_theta{1};
    for i=1:p
        A_theta_bar=A_theta_bar+A_theta{i+1}*theta_bar(i);
    end
    
    % Constraints for every predicted time step
    for k=1:N
        % Dynamics constraints
        
        constraints_MPSC=[constraints_MPSC,x_bar(:,k+1)==(A_theta_bar+B_theta*K)*x_bar(:,k)+B_theta*v(:,k)];
        
        % Tube inclusion constraint
        for i=1:size(H_f,1)
            
            for l=1:size(e_tilde,1)
                matrix_temp=[];
                for j=1:n_MSD-1
                    matrix_temp=[matrix_temp, A_theta{j+1}*x_bar(:,k)];
                end
                constraints_MPSC=[constraints_MPSC,s(:,k+1)>=rho_theta_bar_0*s(:,k)+w_bar+eta*(H_f(i,:)*matrix_temp*e_tilde(l,:)')];
            end
        end
        % State and input constraints
        for j=1:size(F,1)
            constraints_MPSC=[constraints_MPSC,F(j,:)*x_bar(:,k)+G(j,:)*u_bar(:,k)+c(j)*s(:,k)<=1];
        end
        constraints_MPSC=[constraints_MPSC, u_bar(:,k)==v(:,k)+K*x_bar(:,k)];
    end
    % Terminal constraint 
    constraints_MPSC=[constraints_MPSC,x_bar(:,N+1)==zeros(n,1), s(:,N+1)<=1,s(:,N+1)>=0];

    opt_MPSC = optimizer(constraints_MPSC,objective_MPSC,sdpsettings('solver','mosek'),{u_L_k,x_bar(:,1),theta_bar,eta},{u_bar(:,1),x_bar,s});
    
end
%% Feasibility study
% Plot feasibility region
% Get new feasibility optimiser feasibility region
% Objective
objective_MPSC=0;

% Initial constraint
constraints_MPSC=[s(1,1)==0];

A_theta_bar=A_theta{1};
for i=1:p
    A_theta_bar=A_theta_bar+A_theta{i+1}*theta_bar_0(i);
end

% Constraints for every predicted time step
for k=1:N
    % Dynamics constraints
    
    constraints_MPSC=[constraints_MPSC,x_bar(:,k+1)==(A_theta_bar+B_theta*K)*x_bar(:,k)+B_theta*v(:,k)];
    
    % Tube inclusion constraint
    for i=1:size(H_f,1)
        for l=1:size(e_tilde,1)
            matrix_temp=[];
            for j=1:n_MSD-1
                matrix_temp=[matrix_temp, A_theta{j+1}*x_bar(:,k)];
            end
            constraints_MPSC=[constraints_MPSC,s(:,k+1)>=rho_theta_bar_0*s(:,k)+w_bar+eta_0*(H_f(i,:)*matrix_temp*e_tilde(l,:)')];
        end
    end
    % State and input constraints
    for j=1:size(F,1)
        constraints_MPSC=[constraints_MPSC,F(j,:)*x_bar(:,k)+G(j,:)*u_bar(:,k)+c(j)*s(:,k)<=1];
    end
    constraints_MPSC=[constraints_MPSC, u_bar(:,k)==v(:,k)+K*x_bar(:,k)];
end
% Terminal constraint 
constraints_MPSC=[constraints_MPSC,x_bar(:,N+1)==zeros(n,1), s(:,N+1)<=1,s(:,N+1)>=0];

opt_MPSC_feasibility_0 = optimizer(constraints_MPSC,objective_MPSC,sdpsettings('solver','mosek'),{x_bar(:,1)},[]);

%% Simulate MSD from some initial condition
theta_bar_k=theta_bar_0;
eta_k=eta_0;
Omega=theta_bar_0+eta_0*HC;

plot_t=[0];
plot_x=[x_init];
plot_theta_hat=[theta_hat_0];
plot_Omega=[Omega.b];

plot_x_bar=[];
plot_s=[];

term_enlargement_vertices=[];

if ~use_PGSD
    u_L=[];
    for i=1:m
        u_L=[u_L; 3.5*sin(2*pi*(0:T_s:time)+(i-1)*pi/m)];
    end
else
    S=[];
    for i=1:n_MSD
        S=[S;zeros(2,(i-1)),ones(2,1),zeros(2,n_MSD-i)];
    end
    H=10;
    alpha=0.05;
    rng(4)
    theta_PGSD=2*rand(n,m)-1;
    theta_PGSD_orig=theta_PGSD;
end

% Random disturbance sequence (same every time the code is ran)
w_k=rand(n,time/T_s)-0.5;
for i=1:n
    w_k(i,:)=w_k(i,:)*(max(W.V(:,i))-min(W.V(:,i)))-(max(W.V(:,i))+min(W.V(:,i)))/2;
end

x=x_init;
Omega_prev=Omega;
PGSD_counter=0;
u_H=[];
x_H=[];
u_diff=[];
average_time=0;
for k=1:time/T_s
    if use_PGSD
        % Save states
       
        x_H=[x_H, x];
       
        % Perform gradient update
        if PGSD_counter==H
            % Compute gradient wrt to theta_t
            Delta_theta_t=0;
            for t=0:H-1
                Delta_u_t=R*u_H(:,t+1);
                % Compute gradient wrt to u_t 
                for t_1=t+1:H
                    Delta_u_t=Delta_u_t+S'*Q*(x_H(:,t_1+1)-x_ref_PGSD(:,k-H+t_1));
                end
                Delta_theta_t=Delta_theta_t+(x_H(:,t+1)-x_ref_PGSD(:,k-H+t))*Delta_u_t';
            end
            theta_PGSD=theta_PGSD-alpha/H*Delta_theta_t;
            PGSD_counter=0;
            x_H=[x];
            u_H=[];
        end
        % Compute the control input
        u_L(:,k)=theta_PGSD'*(x-x_ref_PGSD(:,k));
        % Save control inputs 
        u_H=[u_H,u_L(:,k)];
        
        % Update the counter
        PGSD_counter=PGSD_counter+1;
    end
    %Compute input
    if use_lorenzen
        
        opt_sol=opt_MPSC(u_L(:,k), x, Omega.A, Omega.b);
    
        u=opt_sol{1}
    else
        tic
        opt_sol=opt_MPSC(u_L(:,k), x, theta_bar_k, eta_k);
        elapsed_time=toc;
        average_time=average_time+elapsed_time;
        u=opt_sol{1};
        u_diff=[u_diff,u-u_L(:,k)];
        term_enlargement_vertices=[term_enlargement_vertices; [opt_sol{2}',opt_sol{3}']];
    end

    % Check if problem is infeasible
    if isnan(u)
        k
        warning('Infeasible!');
        u = K*x;
    end
    
    % Propagate state according to true dynamics with additive disturbance
    x_k1=A*x+B*u+w_k(:,k);
    
    % Collect data for plotting
    plot_t=[plot_t,k*T_s];
    plot_x=[plot_x,x_k1];
    
    
    % Set Membership Update
    % Compute non-falsified parameter set Delta_i
    temp_matrix=[];
    for i=1:p
        temp_matrix=[temp_matrix, A_theta{i+1}*x];
    end
    H_Delta_i=-H_w*temp_matrix;
    h_Delta_i=h_w+H_w*(A_theta{1}*x+B_theta*u-x_k1);
    Delta_i=Polyhedron(H_Delta_i,h_Delta_i);
    
    % Intersection of Theta and Delta
    Intersection=Polyhedron([Omega.A;H_Delta_i],[Omega.b;h_Delta_i]);
    
    
    if (~isEmptySet(Intersection))
        if use_lorenzen
            % SM Update LP
            h=sdpvar(size(h_theta,1),1);
            Lambda=sdpvar(size(h_theta,1),size(h_theta,1)+size(h_w,1),'full');
            objective_SM=ones(1,size(h_theta,1))*h;
            constraints_SM=[Lambda>=0,Lambda*[Omega.A;H_Delta_i]==Omega.A, Lambda*[Omega.b;h_Delta_i]<=h];
            optimize(constraints_SM,objective_SM,sdpsettings('verbose',0));
            
            Omega=Polyhedron(Omega.A,value(h));
        else
            Omega_HB=Intersection.outerApprox();
            
            theta_bar_k1=zeros(p,1);
            eta_k1=zeros(p,1);
            for i=1:p
                theta_bar_k1(i)=(max(Omega_HB.V(:,i))+min(Omega_HB.V(:,i)))/2;
                eta_k1(i)=(max(Omega_HB.V(:,i))-min(Omega_HB.V(:,i)));
            end
            eta_k1=max(eta_k1);
            Projection_Poly=theta_bar_k+(eta_k-eta_k1)*HC;
            eta_k=eta_k1;
            theta_bar_k=Projection_Poly.project(theta_bar_k1).x;
            Omega=theta_bar_k+eta_k*HC;
        end
    else
        warning('Empty Intersection')
    end
    plot_Omega=[plot_Omega, Omega.b];
    
    % Update measurement x
    x=x_k1;
end
u_diff=[u_diff,zeros(m,1)];
average_time=average_time/(time/T_s)

%% Construct nominal matrices for linear MPSC
if compare_linearMPSC
A_lin=A_theta{1};
for i=1:n_MSD-1
     A_lin=A_lin+A_theta{i+1}*((k_min+k_max)/2);%+A_theta{i+n_MSD}*d_true(i);
end
end
Omega_lin=theta_bar_0+eta_0*HC;

% Disturbance set for nominal MPSC
if compare_linearMPSC
    % At every vertex of X, check difference between (A(theta^i)-A_nom)*x^j for 
    % every vertex of Theta
    vertices_W=W.V;
    for i=1:size(X.V,1)
        for j=1:size(Omega_lin.V,1)
            A_err=A_theta{1};
            for k=1:n_MSD-1
                A_err=A_err+A_theta{k+1}*(Omega_lin.V(j,k));%+A_theta{i+n_MSD}*d_true(i);
            end
            A_err=A_err-A_lin;
            vertices_W=[vertices_W;(A_err*X.V(i,:)')'];
        end
    end
    W_lin=Polyhedron(vertices_W).minHRep()
end

%% Compute ellipsoidal RPI set
P_mRPI=sdpvar(n,n);
volume_min=1e4;
lambda=0.7;
    
objective_mRPI=-logdet(P_mRPI);
constraints_mRPI=[];

% Constraint satisfaction
F=[X.A./X.b;zeros(size(U.A,1),n)];
G=[zeros(size(X.A,1),m);U.A./U.b];

for i=1:size(X.A,1)
    constraints_mRPI=[constraints_mRPI, [X.b(i)^2, X.A(i,:); X.A(i,:)',P_mRPI]>=0];
end
for i=1:size(U.A,1)
    constraints_mRPI=[constraints_mRPI, [U.b(i)^2, U.A(i,:)*K; K'*U.A(i,:)',P_mRPI]>=0];
end
for i=1:size(W_lin.V,1)
    constraints_mRPI=[constraints_mRPI, [lambda*P_mRPI-(A_lin+B*K)'*P_mRPI*(A_lin+B*K), -(A_lin+B*K)'*P_mRPI*W_lin.V(i,:)'; -W_lin.V(i,:)*P_mRPI*(A_lin+B*K), -W_lin.V(i,:)*P_mRPI*W_lin.V(i,:)'+1-lambda]>=0];
end
info=optimize(constraints_mRPI, objective_mRPI, sdpsettings('solver', 'mosek', 'verbose', 0))


P_mRPI=value(P_mRPI);

%% Compute support function for minkowski difference
for i=1:size(X.A,1)
 x_set=sdpvar(n,1);
optimize([x_set'*P_mRPI*x_set<=1],-X.A(i,:)*x_set,sdpsettings('solver','mosek','verbose',0));
h_W_X(i)=X.A(i,:)*value(x_set);
end
for i=1:size(U.A,1)
    x_set=sdpvar(n,1);
    optimize([x_set'*P_mRPI*x_set<=1],-U.A(i,:)*K*x_set,sdpsettings('solver','mosek','verbose',0));
    h_W_U(i)=X.A(i,:)*value(x_set);
end
X_f=Polyhedron(H_f,ones(size(H_f,1),1));
for i=1:size(X_f.A,1)
 x_set=sdpvar(n,1);
optimize([x_set'*P_mRPI*x_set<=1],-X_f.A(i,:)*x_set,sdpsettings('solver','mosek','verbose',0));
h_W_X_f(i)=X_f.A(i,:)*value(x_set);
end

%% Set up optimisation for nominal MPSC
% Applied input
u_i=sdpvar(m,N);

% Current learning-based input
u_L_current=sdpvar(m,1);

% Planned safe trajectory
x_i=sdpvar(n,N+1);
x_true=sdpvar(n,1);

% Minimise the norm
objective_MPSC_lin=norm(u_L_current(:,1)-(u_i(:,1)+K*(x_true-x_i(:,1))),2)^2;

constraints_MPSC_lin=[(x_true-x_i(:,1))'*P_mRPI*(x_true-x_i(:,1))<=1];
% For all planned time steps
for i=1:N
    % State propagation
    constraints_MPSC_lin=[constraints_MPSC_lin,x_i(:,i+1)==A_lin*x_i(:,i)+B*u_i(:,i)];
    % State constraints
    constraints_MPSC_lin=[constraints_MPSC_lin,X.A*x_i(:,i)<=X.b-h_W_X'];
    % Input constraints
    constraints_MPSC_lin=[constraints_MPSC_lin, U.A*u_i(:,i)<=U.b-h_W_U'];
end

% Terminal constraint
constraints_MPSC_lin=[constraints_MPSC_lin, x_i(:,N+1)==zeros(n,1)];%<=X_f.b-h_W_X_f'];

% Construct optimiser with input: proposed learning-based control input and
% current state
opt_MPSC_lin=optimizer(constraints_MPSC_lin, objective_MPSC_lin, [], {u_L_current,x_true}, {u_i(:,1)+K*(x_true-x_i(:,1))});

opt_MPSC_lin_feasibility=optimizer(constraints_MPSC_lin, 0, [], {x_true}, []);

%% Simulate MSD for linear MPSC
plot_t_lin=[0];
plot_x_lin=[x_init];


if ~use_PGSD
    u_L=[];
    for i=1:m
        u_L=[u_L; 3.5*sin(2*pi*(0:T_s:time)+(i-1)*pi/m)];
    end
else
    S=[];
    for i=1:n_MSD
        S=[S;zeros(2,(i-1)),ones(2,1),zeros(2,n_MSD-i)];
    end
    H=10;
    alpha=0.05;
    theta_PGSD=theta_PGSD_orig;
end

x=x_init;

PGSD_counter=0;
u_H=[];
x_H=[];
u_diff_lin=[];
for k=1:time/T_s
    if use_PGSD
        % Save states
       
        x_H=[x_H, x];
       
        % Perform gradient update
        if PGSD_counter==H
            % Compute gradient wrt to theta_t
            Delta_theta_t=0;
            for t=0:H-1
                Delta_u_t=R*u_H(:,t+1);
                % Compute gradient wrt to u_t 
                for t_1=t+1:H
                    Delta_u_t=Delta_u_t+S'*Q*(x_H(:,t_1+1)-x_ref_PGSD(:,k-H+t_1));
                end
                Delta_theta_t=Delta_theta_t+(x_H(:,t+1)-x_ref_PGSD(:,k-H+t))*Delta_u_t';
            end
            theta_PGSD=theta_PGSD-alpha/H*Delta_theta_t;
            PGSD_counter=0;
            x_H=[x];
            u_H=[];
        end
        % Compute the control input
        u_L(:,k)=theta_PGSD'*(x-x_ref_PGSD(:,k));
        % Save control inputs 
        u_H=[u_H,u_L(:,k)];
        
        % Update the counter
        PGSD_counter=PGSD_counter+1;
    end
    %Compute input
        opt_sol=opt_MPSC_lin(u_L(:,k), x);
        
        u=opt_sol;
        u_diff_lin=[u_diff_lin,u-u_L(:,k)];

    
    %u=u_L(:,k);
    % Check if problem is infeasible
    if isnan(u)
        k
        warning('Infeasible!');
        u = K*x;
    end
    
    % Propagate state according to true dynamics with additive disturbance
    x_k1=A*x+B*u+w_k(:,k);
    
    % Collect data for plotting
    plot_t_lin=[plot_t_lin,k*T_s];
    plot_x_lin=[plot_x_lin,x_k1];
    
    
    % Update measurement x
    x=x_k1;
end
u_diff_lin=[u_diff_lin,zeros(m,1)];



%% End feasibility study
% Get new feasibility optimiser feasibility region
% Objective
objective_MPSC=0;

% Initial constraint
constraints_MPSC=[s(1,1)==0];

A_theta_bar=A_theta{1};
for i=1:p
    A_theta_bar=A_theta_bar+A_theta{i+1}*theta_bar_k(i);
end

% Constraints for every predicted time step
for k=1:N
    % Dynamics constraints
    
    constraints_MPSC=[constraints_MPSC,x_bar(:,k+1)==(A_theta_bar+B_theta*K)*x_bar(:,k)+B_theta*v(:,k)];
    
    % Tube inclusion constraint
    for i=1:size(H_f,1)
        for l=1:size(e_tilde,1)
            matrix_temp=[];
            for j=1:n_MSD-1
                matrix_temp=[matrix_temp, A_theta{j+1}*x_bar(:,k)];
            end
            constraints_MPSC=[constraints_MPSC,s(:,k+1)>=rho_theta_bar_0*s(:,k)+w_bar+eta_k*(H_f(i,:)*matrix_temp*e_tilde(l,:)')];
        end
    end
    % State and input constraints
    for j=1:size(F,1)
        constraints_MPSC=[constraints_MPSC,F(j,:)*x_bar(:,k)+G(j,:)*u_bar(:,k)+c(j)*s(:,k)<=1];
    end
    constraints_MPSC=[constraints_MPSC, u_bar(:,k)==v(:,k)+K*x_bar(:,k)];
end
% Terminal constraint 
constraints_MPSC=[constraints_MPSC,x_bar(:,N+1)==zeros(n,1), s(:,N+1)<=1,s(:,N+1)>=0];

opt_MPSC_feasibility_k = optimizer(constraints_MPSC,objective_MPSC,sdpsettings('solver','mosek'),{x_bar(:,1)},[]);



%% Terminal Set Enlargement Feasibility
objective_MPSC=0;

% Initial constraint
constraints_MPSC=[s(1,1)==0];

A_theta_bar=A_theta{1};
for i=1:p
    A_theta_bar=A_theta_bar+A_theta{i+1}*theta_bar_k(i);
end

% Constraints for every predicted time step
for k=1:N
    % Dynamics constraints
    
    constraints_MPSC=[constraints_MPSC,x_bar(:,k+1)==(A_theta_bar+B_theta*K)*x_bar(:,k)+B_theta*v(:,k)];
    
    % Tube inclusion constraint
    for i=1:size(H_f,1)
        for l=1:size(e_tilde,1)
            matrix_temp=[];
            for j=1:n_MSD-1
                matrix_temp=[matrix_temp, A_theta{j+1}*x_bar(:,k)];
            end
            constraints_MPSC=[constraints_MPSC,s(:,k+1)>=rho_theta_bar_0*s(:,k)+w_bar+eta_k*(H_f(i,:)*matrix_temp*e_tilde(l,:)')];
        end
    end
    % State and input constraints
    for j=1:size(F,1)
        constraints_MPSC=[constraints_MPSC,F(j,:)*x_bar(:,k)+G(j,:)*u_bar(:,k)+c(j)*s(:,k)<=1];
    end
    constraints_MPSC=[constraints_MPSC, u_bar(:,k)==v(:,k)+K*x_bar(:,k)];
end
% Terminal constraint 
vertices_X=[[zeros(n,1);1],term_enlargement_vertices'];%[[F_0.V,zeros(size(F_0.V,1),1);zeros(1,n),1]',term_enlargement_vertices'];
n_vertices_X=size(vertices_X,2);
lambda_X=sdpvar(n_vertices_X,1,'full');
constraints_MPSC=[constraints_MPSC, sum(lambda_X)==1,lambda_X>=0];
sum_X=zeros(n+1,1);
for i=1:n_vertices_X
    sum_X=sum_X+lambda_X(i)*vertices_X(:,i);
end
constraints_MPSC=[constraints_MPSC, [x_bar(:,N+1);s(:,N+1)]==sum_X];

opt_MPSC_terminal_enlargement = optimizer(constraints_MPSC,objective_MPSC,sdpsettings('solver','mosek'),{x_bar(:,1)},[]);

%%
test=rand(n,1);
tic
sol=opt_MPSC_feasibility_0(test);
toc
tic
sol=opt_MPSC_terminal_enlargement(test);
toc
%%
tic
volume_init=get_feasibility_set_MC(n_MC, n_MSD, X, m, theta_bar_0, eta_0, opt_MPSC_feasibility_0, feasibility_grid_size, true);
toc
tic
volume_end=get_feasibility_set_MC(n_MC, n_MSD, X, m, theta_bar_k, eta_k, opt_MPSC_feasibility_k, feasibility_grid_size, true)
toc

tic
volume_term=get_feasibility_set_MC(n_MC, n_MSD, X, m, theta_bar_k, eta_k, opt_MPSC_terminal_enlargement, feasibility_grid_size, true)
toc

tic
volume_lin=get_feasibility_set_MC(n_MC, n_MSD, X, m, theta_bar_k, eta_k, opt_MPSC_lin_feasibility, feasibility_grid_size, true)
toc

disp('Total Admissible State Space Volume:')
(max(X.V(:,1))-min(X.V(:,1)))^n

%%
if show_figures
    tic
    S_feas{1}=get_feasibility_set(n_MSD, X, m, theta_bar_0, eta_0, opt_MPSC_feasibility_0, feasibility_grid_size, true);
    S_feas{k+1}=get_feasibility_set(n_MSD, X, m, theta_bar_k, eta_k, opt_MPSC_feasibility_k, feasibility_grid_size, true);
    S_feas{k+2}=get_feasibility_set(n_MSD, X, m, theta_bar_k, eta_k, opt_MPSC_terminal_enlargement, feasibility_grid_size, true);
    toc
    show_plots
end

clear objective_mRPI constraints_mRPI x_set x_i x_true u_i u_L_current objective_MPSC_lin constraints_MPSC_lin lambda_X sum_X A_theta_bar eta matrix_temp theta_bar u_bar x_bar constraints_MPSC constraints_PK constraints_SM Delta h h_theta_k H_theta_k Lambda objective_MPSC objective_PK objective_SM options s u_hat u_L_k v x_hat X_PK Y_PK z_X