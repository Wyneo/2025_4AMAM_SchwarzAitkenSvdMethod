function [res_bord, res_mod, list_residu] = Schwarz_Aitken_3c(model1, model2, model3, y0, nb_iter_schwarz, eps, max_iter)
    i=0; 
    y_prec_1=y0{1}-10*eps; %Pas fou mais la flemme de faire mieux
    y_prec_2=y0{2}-10*eps;
    y_prec_3=y0{3}-10*eps;
    list_residu={zeros(1,max_iter-1),zeros(1,max_iter-1),zeros(1,max_iter-1)}; %Résidu pour la première intersection (en haut, sens direct)

    while(i<max_iter && (norm(y_prec_1-y0{1})>eps || norm(y_prec_2-y0{2})>eps || norm(y_prec_3-y0{3})>eps)) %pas sûr du "ou"
        i=i+1;
        [mat_all_iter,cell_all_iter_bord,res_mod_c1,res_mod_c2,res_mod_c3]=iter_solve_3c(model1,model2,model3,nb_iter_schwarz,y0{1}); 
        size(cell_all_iter_bord{1})
        res_bord_c1=acc_aitkenSVD_3c(cell_all_iter_bord{1}); 
        res_bord_c2=acc_aitkenSVD_3c(cell_all_iter_bord{2});
        res_bord_c3=acc_aitkenSVD_3c(cell_all_iter_bord{3});
        y_prec1=y0{1};
        y0{1}=res_bord_c1;
        y_prec2=y0{2};
        y0{2}=res_bord_c2;
        y_prec3=y0{3};
        y0{3}=res_bord_c3;
        list_residu{1}(i)=norm(y_prec1-y0{1});
        list_residu{2}(i)=norm(y_prec2-y0{2});
        list_residu{3}(i)=norm(y_prec3-y0{3});
    end
    list_residu{1}=list_residu{1}(1:i);
    list_residu{2}=list_residu{2}(1:i);
    list_residu{3}=list_residu{3}(1:i);
    res_bord={res_bord_c1, res_bord_c2, res_bord_c3};
    res_mod={res_mod_c1, res_mod_c2, res_mod_c3};
end