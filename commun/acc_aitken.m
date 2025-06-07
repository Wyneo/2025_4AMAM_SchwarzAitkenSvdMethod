function y0=acc_aitken(u)
    yter_now = u(:,2:end);
    yter_prev = u(:,1:end-1);
    err=yter_now - yter_prev;
    err_now = err(:,2:end);
    err_prev = err(:,1:end-1);
    esti_P = err_now * pinv(err_prev);
    y0 = inv((eye(size(esti_P))-esti_P))*(yter_now(:,end)-esti_P*yter_prev(:,end));
end