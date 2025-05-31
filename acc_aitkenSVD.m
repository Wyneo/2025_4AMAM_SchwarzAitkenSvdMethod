function res=acc_aitkenSVD(u)
    [U,S,V]=svd(u);
    u_err_it=u(:,2:end)-u(:,1:end-1);
    u_err_it_now=u_err_it(:,2:end);
    u_err_it_prev=u_err_it(:,1:end-1);
    influence_trace_now=U'*u_err_it_now; 
    influence_trace_prev=U'*u_err_it_prev;
    esti_P=influence_trace_now*pinv(influence_trace_prev); %quasiment matrice nulle ? 
    res=(eye(size(esti_P,1))-esti_P)\(u(:,end)-esti_P*u(:,end-1));
end