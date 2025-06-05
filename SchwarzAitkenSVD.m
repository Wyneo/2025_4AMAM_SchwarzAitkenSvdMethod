function [res_bord, res_mod, list_residu] = SchwarzAitkenSVD(model1, model2, y0, nb_iter_schwarz, eps, max_iter)
    i=0; 
    y_prec=y0-10*eps; 
    list_residu=[10*eps];
    while(i<max_iter && list_residu(end)>eps)
        i=i+1;
        [mat_all_iter,cell_all_iter_bord,res_mod_gauche,res_mod_droit,res_schwarz]=iter_solve(model1,model2,nb_iter_schwarz,y0);
        list_residu=[list_residu; res_schwarz];
        res_bord_droit=acc_aitkenSVD(cell_all_iter_bord{1}); 
        res_bord_gauche=acc_aitkenSVD(cell_all_iter_bord{2});
        y_prec=y0;
        y0=res_bord_droit;
    end
    list_residu=list_residu(find(list_residu>0));
    list_residu=list_residu(2:end); % Enlever la première valeur
    res_bord={res_bord_droit, res_bord_gauche};
    res_mod={res_mod_gauche, res_mod_droit};
end