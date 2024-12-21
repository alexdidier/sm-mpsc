function [A,B] = get_true_Dynamics(n, T_s, k, d)
%GET_TRUE_DYNAMICS 
if any(T_s*d>=0.5) 
    error('Reduce damping!');
end
if n==1
    A=eye(2)+T_s*[0, 1;-k(1), -d(1)];
    B=T_s*[0;1];
else
    A=[T_s*[0, 1, 0, 0; -k(1), -d(1), k(1), d(1)], zeros(2, (n-2)*2)];
    B=T_s*[[0;1], zeros(2,(n-1))];
    for i=2:n-1
        A=[A; zeros(2, (i-2)*2), T_s*[0,0,0,1,0,0; k(i-1), d(i-1), -(k(i-1)+k(i)), -(d(i-1)+d(i)), k(i), d(i)], zeros(2,((n-1)-i)*2)];
        B=[B;zeros(2,i-1),T_s*[0;1],zeros(2, n-i)];
    end
    A=[A; zeros(2,(n-2)*2),T_s*[0, 0, 0, 1; k(n-1), d(n-1), -k(n-1), -d(n-1)]];
    A=eye(2*n)+A;
    B=[B;zeros(2,n-1), T_s*[0;1]];
end

