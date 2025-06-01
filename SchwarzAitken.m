function [res_bord, res_mod, list_residu] = Schwarz_Aitken(model1, model2, y0, nb_iter_schwarz, eps, max_iter)
    i=0; 
    y_prec=y0-10*eps; %Pas fou mais la flemme de faire mieux
    list_residu=zeros(1,max_iter-1);
    while(i<max_iter && norm(y_prec-y0)>eps)
        i=i+1;
        [mat_all_iter,cell_all_iter_bord,res_mod_gauche,res_mod_droit]=iter_solve(model1,model2,nb_iter_schwarz,y0); 
        res_bord_droit=acc_aitkenSVD(cell_all_iter_bord{1}); 
        res_bord_gauche=acc_aitkenSVD(cell_all_iter_bord{2});
        y_prec=y0;
        y0=res_bord_droit;
        list_residu(i)=norm(y_prec-y0);
    end
    list_residu=list_residu(1:i);
    res_bord={res_bord_droit, res_bord_gauche};
    res_mod={res_mod_gauche, res_mod_droit};
end