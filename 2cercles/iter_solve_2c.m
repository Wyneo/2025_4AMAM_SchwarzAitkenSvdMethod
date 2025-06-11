function [all_iter,all_iter_bord, results1, results2, res_schwarz, bool_convergence]=iter_solve_2c(model1,model2,itemax,y0,eps) %Itérations sur les deux sous domaines
    %model1,model2 : les sous domaines de la géométrie
    %itemax : nombre d'itérations voulues 
    %y0 : condition initiale sur le bord du premier sous domaine

    cd_1=findNodes(model1.Mesh,"region","Edge",[2,3]); %arc de cercle à droite du premier sous domaine
    cg_2=findNodes(model2.Mesh,"region","Edge",[1,5]); %arc de cercle à gauche du deuxième sous domaine

    all_iter_cd=[]; 
    all_iter_cg=[];
    all_iter_bord_cd=[y0{1}]; 
    all_iter_bord_cg=[y0{2}];
    res_schwarz=[];
    res_schwarz_1=[];
    res_schwarz_2=[];

    %Initialisation
    % Première itération pour "intégrer le y0" dans la solution
    applyBoundaryCondition(model1,"dirichlet","Edge",[1,4:7],"u",0);
    applyBoundaryCondition(model1,"dirichlet","Edge",[2,3],"u",0);

    u_temp=solvepde(model1).NodalSolution;
    u_temp(cd_1)=y0{1};

    [p1,e1,t1]=meshToPet(model1.Mesh);
    F1 = pdeInterpolant(p1,t1,u_temp);
    cl_cd_1=@(region,state) evaluate(F1,region.x,region.y);

    i=1; % Compteur d'itérations
    bool_convergence=false;
    u1=u_temp; 
    u2=u_temp; 
    for i=1:itemax
        %Premier sous domaine
        stock_u1_res=u1(cd_1);
        applyBoundaryCondition(model1,"dirichlet","Edge",[1,4:7],"u",0);
        applyBoundaryCondition(model1,"dirichlet","Edge",[2,3],"u",cl_cd_1);

        results1=solvepde(model1);
        u1=results1.NodalSolution;

        %res_schwarz_mod1=norm(u1(cd_1)-y_sol{1},Inf);
        res_schwarz_1=norm(u1(cd_1)-stock_u1_res,Inf);

        %Récupération de la solution sur le bord
        [p1,e1,t1]=meshToPet(model1.Mesh);
        F1 = pdeInterpolant(p1,t1,u1);
        cl_cg_2=@(region,state) evaluate(F1,region.x,region.y); %condition initiale

        %Deuxième sous domaine
        stock_u2_res=u2(cg_2);
        applyBoundaryCondition(model2,"dirichlet","Edge",[2:4,6:7],"u",0);
        applyBoundaryCondition(model2,"dirichlet","Edge",[1,5],"u",cl_cg_2);

        results2=solvepde(model2);
        u2=results2.NodalSolution;

        res_schwarz_2=norm(u2(cg_2)-stock_u2_res,Inf);

        [p2,e2,t2]=meshToPet(model2.Mesh);
        F2 = pdeInterpolant(p2,t2,u2);
        cl_cd_1=@(region,state) evaluate(F2,region.x,region.y); %condition initiale

        %Stockage des itérations
        all_iter_cd=[all_iter_cd,u1];
        all_iter_cg=[all_iter_cg,u2];
        all_iter_bord_cd=[all_iter_bord_cd,u1(cd_1)]; % bord droit
        all_iter_bord_cg=[all_iter_bord_cg,u2(cg_2)]; % bord gauche
        res_schwarz=[res_schwarz; max(res_schwarz_1,res_schwarz_2)]; % Résidu de Schwarz
        %disp(res_schwarz(end)); 
        i=i+1;
    end
    res_schwarz=res_schwarz(2:end); % Enlever la première valeur 
    bool_convergence=(res_schwarz(end)<eps);
    all_iter_bord={all_iter_bord_cd,all_iter_bord_cg}; % cell pour pouvoir accéder aux bords séparément
    all_iter={all_iter_cd,all_iter_cg}; % cell pour pouvoir accéder aux itérations séparément
end