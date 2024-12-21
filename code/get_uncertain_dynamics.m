function [A_theta] = get_uncertain_dynamics(n, T_s, k_only, d)
%GET_UNCERTAIN_DYNAMICS 
A_theta{1}=[T_s*[0, 1, 0, 0; 0, 0, 0, 0], zeros(2, (n-2)*2)];
for i=2:n-1
    A_theta{1}=[A_theta{1}; zeros(2, (i-2)*2), T_s*[0,0,0,1,0,0; 0, 0, 0, 0,0,0], zeros(2,((n-1)-i)*2)];
end
A_theta{1}=[A_theta{1}; zeros(2,(n-2)*2),T_s*[0, 0, 0, 1; 0, 0, 0, 0]];
if k_only
    for i=n:2*(n-1)
        j=i-n+1;
        A_theta{1}=A_theta{1}+d(j)*[zeros(2*(j-1),2*n); zeros(4, (j-1)*2), T_s*[zeros(1,4); 0,-1,0,1; zeros(1,4);0,1,0,-1], zeros(4, (n-1-j)*2); zeros(2*(n-1-j), 2*n)];
    end
end
A_theta{1}=eye(2*n)+A_theta{1};

for i=1:n-1
    A_theta{i+1}=[zeros(2*(i-1),2*n); zeros(4, (i-1)*2), T_s*[zeros(1,4); -1,0,1,0; zeros(1,4);1,0,-1,0], zeros(4, (n-1-i)*2); zeros(2*(n-1-i), 2*n)];
end
if ~k_only
    for i=n:2*(n-1)
        j=i-n+1;
        A_theta{i+1}=[zeros(2*(j-1),2*n); zeros(4, (j-1)*2), T_s*[zeros(1,4); 0,-1,0,1; zeros(1,4);0,1,0,-1], zeros(4, (n-1-j)*2); zeros(2*(n-1-j), 2*n)];
    end
end

end

