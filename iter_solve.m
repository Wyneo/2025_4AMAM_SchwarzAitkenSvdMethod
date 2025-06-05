function [all_iter,all_iter_bord, results1, results2, res_schwarz]=iter_solve(model1,model2,itemax,y0) %Itérations sur les deux sous domaines
    %model1,model2 : les sous domaines de la géométrie
    %itemax : nombre d'itérations voulues 

    cd_1=findNodes(model1.Mesh,"region","Edge",3); %arc de cercle à droite du premier sous domaine
    cg_2=findNodes(model2.Mesh,"region","Edge",5); %arc de cercle à gauche du deuxième sous domaine

    cl_cd_1 = y0;

    all_iter_cd=[]; 
    all_iter_cg=[];
    all_iter_bord_cd=[]; 
    all_iter_bord_cg=[];
    res_schwarz=[];

    %Initialisation
    % Première itération pour "intégrer le y0" dans la solution
    applyBoundaryCondition(model1,"dirichlet","Edge",[1,2,4,5],"u",0);
    applyBoundaryCondition(model1,"dirichlet","Edge",3,"u",0);

    u_temp=solvepde(model1).NodalSolution;
    u_temp(cd_1)=y0;

    [p1,e1,t1]=meshToPet(model1.Mesh);
    F1 = pdeInterpolant(p1,t1,u_temp);
    cl_cd_temp=@(region,state) evaluate(F1,region.x,region.y);

    %Premier sous domaine
    applyBoundaryCondition(model1,"dirichlet","Edge",[1,2,4,5],"u",0);
    applyBoundaryCondition(model1,"dirichlet","Edge",3,"u",cl_cd_temp);

    results1=solvepde(model1);
    u1=results1.NodalSolution;

    %Récupération de la solution sur le bord
    [p1,e1,t1]=meshToPet(model1.Mesh);
    F1 = pdeInterpolant(p1,t1,u1);
    cl_cg_2=@(region,state) evaluate(F1,region.x,region.y); %condition initiale

    %Deuxième sous domaine
    applyBoundaryCondition(model2,"dirichlet","Edge",[1,2,3,4],"u",0);
    applyBoundaryCondition(model2,"dirichlet","Edge",5,"u",cl_cg_2);

    results2=solvepde(model2);
    u2=results2.NodalSolution;

    [p2,e2,t2]=meshToPet(model2.Mesh);
    F2 = pdeInterpolant(p2,t2,u2);
    cl_cd_1=@(region,state) evaluate(F2,region.x,region.y); %condition initiale

    %Stockage des itérations
    all_iter_cd=[all_iter_cd,u1];
    all_iter_cg=[all_iter_cg,u2];
    all_iter_bord_cd=[all_iter_bord_cd,u1(cd_1)]; % bord droit
    all_iter_bord_cg=[all_iter_bord_cg,u2(cg_2)]; % bord gauche

    for i=2:itemax
        %Premier sous domaine
        stock_u_res=u1(cd_1);
        applyBoundaryCondition(model1,"dirichlet","Edge",[1,2,4,5],"u",0);
        applyBoundaryCondition(model1,"dirichlet","Edge",3,"u",cl_cd_1);

        results1=solvepde(model1);
        u1=results1.NodalSolution;

        res_schwarz=[res_schwarz;norm(u1(cd_1)-stock_u_res)];

        %Récupération de la solution sur le bord
        [p1,e1,t1]=meshToPet(model1.Mesh);
        F1 = pdeInterpolant(p1,t1,u1);
        cl_cg_2=@(region,state) evaluate(F1,region.x,region.y); %condition initiale

        %Deuxième sous domaine
        applyBoundaryCondition(model2,"dirichlet","Edge",[1,2,3,4],"u",0);
        applyBoundaryCondition(model2,"dirichlet","Edge",5,"u",cl_cg_2);

        results2=solvepde(model2);
        u2=results2.NodalSolution;

        [p2,e2,t2]=meshToPet(model2.Mesh);
        F2 = pdeInterpolant(p2,t2,u2);
        cl_cd_1=@(region,state) evaluate(F2,region.x,region.y); %condition initiale

        %Stockage des itérations
        all_iter_cd=[all_iter_cd,u1];
        all_iter_cg=[all_iter_cg,u2];
        all_iter_bord_cd=[all_iter_bord_cd,u1(cd_1)]; % bord droit
        all_iter_bord_cg=[all_iter_bord_cg,u2(cg_2)]; % bord gauche
    end

    all_iter_bord={all_iter_bord_cd,all_iter_bord_cg}; % cell pour pouvoir accéder aux bords séparément
    all_iter={all_iter_cd,all_iter_cg}; % cell pour pouvoir accéder aux itérations séparément
end