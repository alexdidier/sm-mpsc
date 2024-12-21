function Feasibility = get_feasibility_set(n_MSD, X, m, theta_bar, eta, opt_MPSC, grid_size, verbose)
%GET_FEASIBILITY_SET 

if nargin<7
    grid_size=0.1;
    verbose=true;
elseif nargin==7
    verbose=true;
end

feasible=[];
if verbose
   prog=0;
   prog_max=length(min(X.V(:,1)):grid_size:max(X.V(:,1)));
   fprintf(1, 'Feasibility Analysis Progress: %3d%%\n', prog) 
end




for x_1=min(X.V(:,1)):grid_size:max(X.V(:,1))
    if verbose
        prog=(x_1-min(X.V(:,1)))/grid_size/prog_max*100;
        fprintf(1, '\b\b\b\b%3.0f%%', prog);
    end
    for x_2=min(X.V(:,3)):grid_size:max(X.V(:,3))
        for x_3=min(X.V(:,5)):grid_size:max(X.V(:,5))
            x=[x_1;2;x_2;-2;x_3;2];
            if n_MSD>3
                x=[x;zeros(2*(n_MSD-3),1)];
            end
            [opt_sol, errorcode]=opt_MPSC(x);
            if errorcode==0
                 feasible=[feasible;x'];
            end
        end
    end
end
if verbose
    fprintf(1, '\b\b\b\b%3.0f%%\n', 100);
end
fprintf(1,'Computing Feasibility Polyhedron:\n');
fprintf(1,'Size of feasible points:\n');
disp(size(feasible))
Feasibility=Polyhedron(feasible(:,[1,3,5]));
Feasibility.minVRep();
end


