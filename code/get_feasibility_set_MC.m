function volume = get_feasibility_set_MC(N, n_MSD, X, m, theta_bar, eta, opt_MPSC, grid_size, verbose)
%GET_FEASIBILITY_SET 
if nargin<7
    grid_size=0.1;
    verbose=true;
elseif nargin==7
    verbose=true;
end


input=zeros(m,1);
feasible=0;

x_1=min(X.V(:,1)):grid_size:max(X.V(:,1));
scaling=max(X.V(:,1))-min(X.V(:,1));
offset=min(X.V(:,1));

parfor i=1:N
    [opt_sol, errorcode]=opt_MPSC(rand(2*n_MSD,1)*scaling+offset);
    if errorcode==0
        feasible=feasible+1;
    end
end

fprintf(1,'Number of feasible points:\n');
disp(feasible);

volume=scaling^(2*n_MSD)/N*feasible;

fprintf(1,'Approximate Polyhedron Volume:');
disp(volume);
end


