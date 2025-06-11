function [res_bord, res_mod, list_residu] = SchwarzAitkenSVD_3c(model1, model2, model3, y0, nb_iter_schwarz, eps, max_iter)
    i=0; 
    y_prec={y0{1}-1,y0{2}-1,y0{3}-1};
    list_residu=[10]; %Résidu pour la première interface seulement (en haut)
    bool_convergence=false; % Pour arrêter si le résidu augmente
    while(i<max_iter && list_residu(end)>eps && not(bool_convergence))
        i=i+1;
        [mat_all_iter,cell_all_iter_bord,res_mod_c1,res_mod_c2,res_mod_c3,res_schwarz,bool_convergence]=iter_solve_3c(model1,model2,model3,nb_iter_schwarz,y0,eps);
        list_residu=[list_residu; res_schwarz];
        res_bord_c1=acc_aitkenSVD(cell_all_iter_bord{1}, 1e-5); 
        res_bord_c2=acc_aitkenSVD(cell_all_iter_bord{2}, 1e-5);
        res_bord_c3=acc_aitkenSVD(cell_all_iter_bord{3}, 1e-5);
        y_prec={y0{1},y0{2},y0{3}};
        y0={res_bord_c1,res_bord_c2,res_bord_c3};
        if(list_residu(end)-list_residu(end-(nb_iter_schwarz-2))>eps)
            disp("Augmentation du résidu, arrêt des itérations : SVD");
            bool_convergence=true; % Dans le cas où le résidu augmente, on arrête
        end
    end
    list_residu=list_residu(2:end); % Enlever la première valeur
    res_bord={res_bord_c1, res_bord_c2, res_bord_c3};
    res_mod={res_mod_c1, res_mod_c2, res_mod_c3};
end