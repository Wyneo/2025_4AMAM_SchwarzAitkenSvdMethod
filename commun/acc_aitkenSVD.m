function y0=acc_aitkenSVD(u)
    [U,S,V]=svd(u);
    epsilon=1e-10;
    diag_S = diag(S)
    ind=find(diag_S>epsilon/diag_S(1));
    n_gamma = length(ind);
    %! Tentative de méthode du coude
    % vec_courbure = [diag_S(1)-5*diag_S(2)+4*diag_S(3)-diag_S(4);diag_S(1:end-2) - 2*diag_S(2:end-1) + diag_S(3:end);diag_S(end)-5*diag_S(end-1)+4*diag_S(end-2)-diag_S(end-3)];
    % [t_max,ind]=max(vec_courbure); %maximum de la courbure
    % disp(ind)
    U=U(:,ind);
    S=S(ind,ind);
    V=V(:,ind);
    y_chap_now=U'*u(:,end-n_gamma:end);
    y_chap_prev=U'*u(:,end-n_gamma-1:end-1);
    e_chap=y_chap_now-y_chap_prev;
    esti_P=e_chap(:,2:end)*pinv(e_chap(:,1:end-1));
    y0 = U*inv((eye(n_gamma)-esti_P))*(y_chap_now(:,end)-esti_P*y_chap_prev(:,end));
end