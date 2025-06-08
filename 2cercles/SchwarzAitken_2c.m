function [res_bord, res_mod, list_residu] = SchwarzAitken_2c(model1, model2, y0, y_sol, nb_iter_schwarz, eps, max_iter)
    i=0; 
    list_residu=[10];
    bool_convergence=false;
    while(i<max_iter && list_residu(end)>eps && not(bool_convergence))
        i=i+1;
        [mat_all_iter,cell_all_iter_bord,res_mod_gauche,res_mod_droit,res_schwarz,bool_convergence]=iter_solve_2c(model1,model2,nb_iter_schwarz,y0,y_sol,eps);
        list_residu=[list_residu; res_schwarz];
        res_bord_droit=acc_aitken(cell_all_iter_bord{1}); 
        res_bord_gauche=acc_aitken(cell_all_iter_bord{2});
        y_prec=y0;
        y0=res_bord_droit;
        if(list_residu(end)-list_residu(end-(nb_iter_schwarz-1))>eps)
            disp("Augmentation du résidu, arrêt des itérations : Simple ");
            bool_convergence=true; % Dans le cas où le résidu augmente, on arrête
        end
    end
    list_residu = list_residu(2:end); % Enlever la première valeur
    res_bord={res_bord_droit, res_bord_gauche};
    res_mod={res_mod_gauche, res_mod_droit};
end