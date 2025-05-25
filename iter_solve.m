function all_iter=iter_solve(model1,model2,itemax) %Itérations sur les deux sous domaines
    %model1,model2 : les sous domaines de la géométrie
    %itemax : nombre d'itérations voulues 

    cd_1=findNodes(model1.Mesh,"region","Edge",3); %arc de cercle à droite du premier sous domaine
    cg_2=findNodes(model2.Mesh,"region","Edge",5); %arc de cercle à gauche du deuxième sous domaine
    cl_cd_1=0; %condition initiale
    all_iter=[];

    for i=1:itemax

        %Premier sous domaine
        applyBoundaryCondition(model1,"dirichlet","Edge",[1,2,4,5],"u",0);
        applyBoundaryCondition(model1,"dirichlet","Edge",3,"u",cl_cd_1);

        results1=solvepde(model1);
        u1=results1.NodalSolution;

        %Récupération de la solution sur le bord
        [p1,e1,t1]=meshToPet(model1.Mesh);
        F1 = pdeInterpolant(p1,t1,u1);
        cl_cg_2=@(region,state) evaluate(F1,region.x,region.y); %condition initiale

        %Deuxième sous domaine
        applyBoundaryCondition(model2,"dirichlet","Edge",[1,2,3,4],"u",0);
        applyBoundaryCondition(model2,"dirichlet","Edge",[5],"u",cl_cg_2);

        results2=solvepde(model2);
        u2=results2.NodalSolution;

        [p2,e2,t2]=meshToPet(model2.Mesh);
        F2 = pdeInterpolant(p2,t2,u2);
        cl_cd_1=@(region,state) evaluate(F2,region.x,region.y); %condition initiale

        %Stockage des itérations
        all_iter=[all_iter [u1;u2]];
    end
end