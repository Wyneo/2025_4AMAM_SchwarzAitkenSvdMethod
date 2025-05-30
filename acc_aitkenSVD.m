function res=acc_aitkenSVD(u,epsilon)
    [U,S,V]=svd(u);
    sig_vec=diag(S);
    bound_inf=epsilon*S(1,1); % see draft with circles
    [temp,ind_sig_vec_min]=min(abs(sig_vec-bound_inf)); % Detect the nearest of index bound_inf. 
    U_rest=U(1:ind_sig_vec_min,1:ind_sig_vec_min);
    S_rest=S(1:ind_sig_vec_min,1:ind_sig_vec_min);
    V_rest=V(1:ind_sig_vec_min,1:ind_sig_vec_min); 
    u_err=u(:,2:end)-u(:,1:end-1);
    estim_p=u_err*U_rest'*pinv(u_err*U_rest'); 
    disp(size(estim_p))
    y_inf=(eye(ind_sig_vec_min,ind_sig_vec_min)-estim_p)*(U(:,ind_sig_vec_min+2)-estim_p*U(:,ind_sig_vec_mi+1))
    res=y_inf;
end