function [k,d] = get_MSD_params(n_MSD, k_min,k_max,d_min,d_max)
%GET_MSD_PARAMS 
    k=zeros(n_MSD-1,1);
    d=zeros(n_MSD-1,1);
    for i=1:n_MSD-1
        k(i)=rand()*(k_max-k_min)+k_min;
        d(i)=rand()*(d_max-d_min)+d_min;
    end
end

