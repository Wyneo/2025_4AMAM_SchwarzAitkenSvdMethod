function res=acc_aitkenSVD(u,epsilon)
    [U,S,V]=svd(acc_aitkenSVD);
    sig_vec=diag(S);
    bound_inf=epsilon*S(1,1); % see draft with circles
    [temp,ind_sig_vec_min]=which.min(abs(sig_vec-bound_inf)); % Detect the nearest of index bound_inf. 
    [U_rest,S_rest,V_rest]=[U(1:ind_sig_vec_min,1:ind_sig_vec_min),
                            S(1:ind_sig_vec_min,1:ind_sig_vec_min),
                            V(1:ind_sig_vec_min,1:ind_sig_vec_min)]; 
    estim_p=S_rest; 
    y_inf=(eye(ind_sig_vec_min,ind_sig_vec_min)-estim_p)*(U(:,ind_sig_vec_min+2)-estim_p*U(:,ind_sig_vec_min))
    res=y_inf;
end