function [res_bord, res_mod, list_residu] = SchwarzAitken_3c(model1, model2, model3, y0, nb_iter_schwarz, eps, max_iter)
    i=0; 
    list_residu=[10];
    bool_convergence=false;
    while(i<max_iter && list_residu(end)>eps && not(bool_convergence))
        i=i+1;
        [mat_all_iter,cell_all_iter_bord,res_mod_c1,res_mod_c2,res_mod_c3,res_schwarz,bool_convergence]=iter_solve_3c(model1,model2,model3,nb_iter_schwarz,y0,eps);
        list_residu=[list_residu; res_schwarz];
        res_bord_c1=acc_aitken(cell_all_iter_bord{1}); 
        res_bord_c2=acc_aitken(cell_all_iter_bord{2});
        res_bord_c3=acc_aitken(cell_all_iter_bord{3});
        y_prec=y0;
        y0=res_bord_c1;
        if(list_residu(end)-list_residu(end-(nb_iter_schwarz-2))>eps)
            disp("Augmentation du résidu, arrêt des itérations : Simple ");
            bool_convergence=true; % Dans le cas où le résidu augmente, on arrête
        end
    end
    list_residu = list_residu(2:end); % Enlever la première valeur
    res_bord={res_bord_c1, res_bord_c2, res_bord_c3};
    res_mod={res_mod_c1, res_mod_c2, res_mod_c3};
end