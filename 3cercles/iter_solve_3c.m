function [all_iter,all_iter_bord, results1, results2, results3, res_schwarz, bool_convergence]=iter_solve_3c(model1,model2,model3,itemax,y0,eps) %Itérations sur les deux sous domaines
    %model1,model2,model3 : les sous domaines de la géométrie
    %itemax : nombre d'itérations voulues 
    %y0 : condition initiale sur le bord du premier sous domaine
    
    c1=findNodes(model1.Mesh,"region","Edge",[2,3,4]);
    c2=findNodes(model2.Mesh,"region","Edge",[6,1,2]);
    c3=findNodes(model3.Mesh,"region","Edge",[3,4,5]);

    all_iter_c1=[];
    all_iter_c2=[];
    all_iter_c3=[];
    all_iter_bord_c1=[];
    all_iter_bord_c2=[];
    all_iter_bord_c3=[];
    res_schwarz=[];
    res_schwarz_1=[];
    res_schwarz_2=[];
    res_schwarz_3=[];

    %Initialisation
    % Première itération pour "intégrer le y0" dans la solution
    applyBoundaryCondition(model1,"dirichlet","Edge",[1,5:10],"u",0);
    applyBoundaryCondition(model1,"dirichlet","Edge",[2,3,4],"u",0);

    u_temp=solvepde(model1).NodalSolution;
    u_temp(c1)=y0;

    [p1,e1,t1]=meshToPet(model1.Mesh);
    F1 = pdeInterpolant(p1,t1,u_temp);
    cl_c1temp=@(region,state) evaluate(F1,region.x,region.y);

    %Premier sous domaine
    applyBoundaryCondition(model1,"dirichlet","Edge",[1,5:10],"u",0);
    applyBoundaryCondition(model1,"dirichlet","Edge",[2,3,4],"u",cl_c1temp);

    results1=solvepde(model1);
    u1=results1.NodalSolution;

    %Récupération de la solution sur le bord
    [p1,e1,t1]=meshToPet(model1.Mesh);
    F1 = pdeInterpolant(p1,t1,u1);
    u_c1=@(region,state) evaluate(F1,region.x,region.y); %résultat sur bord du premier sous domaine

    %Deuxième sous domaine
    applyBoundaryCondition(model2,"dirichlet","Edge",[2:5,7:10],"u",0);
    applyBoundaryCondition(model2,"dirichlet","Edge",[1,6],"u",u_c1);

    results2=solvepde(model2);
    u2=results2.NodalSolution;

    %Récupération de la solution sur le bord
    [p2,e2,t2]=meshToPet(model2.Mesh);
    F2 = pdeInterpolant(p2,t2,u2);
    u_c2=@(region,state) evaluate(F2,region.x,region.y); %résultat sur bord du deuxième sous domaine

    %Troisième sous domaine
    applyBoundaryCondition(model3,"dirichlet","Edge",[1,2,6:10],"u",0);
    applyBoundaryCondition(model3,"dirichlet","Edge",[3,4],"u",u_c2);
    applyBoundaryCondition(model3,"dirichlet","Edge",[5,4],"u",u_c1);

    results3=solvepde(model3);
    u3=results3.NodalSolution;

    %Récupération de la solution sur le bord
    [p3,e3,t3]=meshToPet(model3.Mesh);
    F3 = pdeInterpolant(p3,t3,u3);
    u_c3=@(region,state) evaluate(F3,region.x,region.y); %résultat sur bord du troisième sous domaine

    %Stockage des itérations
    all_iter_c1=[all_iter_c1,u1];
    all_iter_c2=[all_iter_c2,u2];
    all_iter_c3=[all_iter_c3,u3];
    all_iter_bord_c1=[all_iter_bord_c1,u1(c1)]; % bord premier sous domaine
    all_iter_bord_c2=[all_iter_bord_c2,u2(c2)]; % bord deuxième sous domaine
    all_iter_bord_c3=[all_iter_bord_c3,u3(c3)]; % bord troisième sous domaine
    
    i=1; % Compteur d'itérations
    bool_convergence=false;
    u1=u_temp; 
    u2=u_temp; 
    for i=2:itemax
        %Premier sous domaine
        stock_u1_res=u1(c1);
        applyBoundaryCondition(model1,"dirichlet","Edge",[1,5:10],"u",0);
        applyBoundaryCondition(model1,"dirichlet","Edge",4,"u",u_c2);
        applyBoundaryCondition(model1,"dirichlet","Edge",[2,3],"u",u_c3);

        results1=solvepde(model1);
        u1=results1.NodalSolution;

        res_schwarz_1=norm(u1(c1)-stock_u1_res,Inf);

        %Récupération de la solution sur le bord
        [p1,e1,t1]=meshToPet(model1.Mesh);
        F1 = pdeInterpolant(p1,t1,u1);
        u_c1=@(region,state) evaluate(F1,region.x,region.y); %résultat sur bord du premier sous domaine

        %Deuxième sous domaine
        stock_u2_res=u2(c2);
        applyBoundaryCondition(model2,"dirichlet","Edge",[3:5,7:10],"u",0);
        applyBoundaryCondition(model2,"dirichlet","Edge",2,"u",u_c3);
        applyBoundaryCondition(model2,"dirichlet","Edge",[1,6],"u",u_c1);

        results2=solvepde(model2);
        u2=results2.NodalSolution;

        %Récupération de la solution sur le bord
        [p2,e2,t2]=meshToPet(model2.Mesh);
        F2 = pdeInterpolant(p2,t2,u2);
        u_c2=@(region,state) evaluate(F2,region.x,region.y); %résultat sur bord du deuxième sous domaine

        res_schwarz_2=norm(u2(c2)-stock_u2_res,Inf);

        %Troisième sous domaine
        stock_u3_res=u3(c3);
        applyBoundaryCondition(model3,"dirichlet","Edge",[1,2,6:10],"u",0);
        applyBoundaryCondition(model3,"dirichlet","Edge",5,"u",u_c1);
        applyBoundaryCondition(model3,"dirichlet","Edge",[3,4],"u",u_c2);
        results3=solvepde(model3);
        u3=results3.NodalSolution;

        %Récupération de la solution sur le bord
        [p3,e3,t3]=meshToPet(model3.Mesh);
        F3 = pdeInterpolant(p3,t3,u3);
        u_c3=@(region,state) evaluate(F3,region.x,region.y); %résultat sur bord du troisième sous domaine

        res_schwarz_3=norm(u3(c3)-stock_u3_res,Inf);

        %Stockage des itérations
        all_iter_c1=[all_iter_c1,u1];
        all_iter_c2=[all_iter_c2,u2];
        all_iter_c3=[all_iter_c3,u3];
        all_iter_bord_c1=[all_iter_bord_c1,u1(c1)]; % bord premier sous domaine
        all_iter_bord_c2=[all_iter_bord_c2,u2(c2)]; % bord deuxième sous domaine
        all_iter_bord_c3=[all_iter_bord_c3,u3(c3)]; % bord troisième sous domaine
        i=i+1;
        res_schwarz=[res_schwarz; max(max(res_schwarz_1,res_schwarz_2),res_schwarz_3)]; % Résidu de Schwarz
    end
    res_schwarz=res_schwarz(2:end); % Enlever la première valeur 
    bool_convergence=(res_schwarz(end)<eps);
    all_iter_bord={all_iter_bord_c1,all_iter_bord_c2,all_iter_bord_c3}; % cell pour pouvoir accéder aux bords séparément
    all_iter={all_iter_c1,all_iter_c2,all_iter_c3}; % cell pour pouvoir accéder aux itérations séparément
end