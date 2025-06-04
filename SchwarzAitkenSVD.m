function [res_bord, res_mod, list_residu] = SchwarzAitkenSVD(model1, model2, y0, nb_iter_schwarz, eps, max_iter)
    i=0; 
    y_prec=y0-10*eps; %Pas fou mais la flemme de faire mieux
    list_residu=zeros(max_iter*(nb_iter_schwarz),1);
    while(i<max_iter && norm(y_prec-y0)>eps)
        i=i+1;
        [mat_all_iter,cell_all_iter_bord,res_mod_gauche,res_mod_droit,res_schwarz]=iter_solve(model1,model2,nb_iter_schwarz,y0);
        list_residu(1+i*nb_iter_schwarz:i*nb_iter_schwarz+nb_iter_schwarz-1)=res_schwarz;
        res_bord_droit=acc_aitkenSVD(cell_all_iter_bord{1}); 
        res_bord_gauche=acc_aitkenSVD(cell_all_iter_bord{2});
        y_prec=y0;
        y0=res_bord_droit;
        list_residu((i+1)*nb_iter_schwarz)=norm(y_prec-y0);
    end
    list_residu=list_residu(find(list_residu>0));
    res_bord={res_bord_droit, res_bord_gauche};
    res_mod={res_mod_gauche, res_mod_droit};
end