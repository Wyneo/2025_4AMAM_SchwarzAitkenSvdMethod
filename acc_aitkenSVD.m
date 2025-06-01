function y_0=acc_aitkenSVD(u)
    [U,S,V]=svd(u);
    epsilon=1e-10;
    ind=find(diag(S)>epsilon);
    n_gamma = length(ind);
    U=U(:,ind);
    S=S(ind,ind);
    V=V(:,ind);
    y_chap_now=U'*u(:,end-n_gamma:end);
    y_chap_prev=U'*u(:,end-n_gamma-1:end-1);
    e_chap=y_chap_now-y_chap_prev;
    esti_P=e_chap(:,2:end)*pinv(e_chap(:,1:end-1));
    y_0 = U*inv((eye(n_gamma)-esti_P))*(y_chap_now(:,end)-esti_P*y_chap_prev(:,end-1));
end