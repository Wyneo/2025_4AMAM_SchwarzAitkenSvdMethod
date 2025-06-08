clear 
close all

addpath("./commun")

R=2;
L=2; 
y=0; 
x=0;

circle2_sol=circle2_solution_manufactured(R,L);

long_EF_max=0.2;
ord_EF="Quadratic";

model1=create1circleintermodel(x,y,L,R,long_EF_max,ord_EF);
model2=create1circleintermodel(x,y,-L,R,long_EF_max,ord_EF);

u_fun=circle2_solution_manufactured(R,L);

% ! Tentative d'appliquer l'opérateur de Laplace Discret, puis d'utiliser la solution pour définir la source

u_sol_circle1=u_fun(model1.Mesh.Nodes(1,:),model1.Mesh.Nodes(2,:));
u_sol_circle2=u_fun(model2.Mesh.Nodes(1,:),model2.Mesh.Nodes(2,:));

figure(1)
pdeplot(model1.Mesh,"XYData",u_sol_circle1)
saveas(gcf,"Sol_Circle1.jpg")
figure(2)
pdeplot(model2.Mesh,"XYData",u_sol_circle2)
saveas(gcf,"Sol_Circle2.jpg")

specifyCoefficients(model1,"m",0,"d",0,"c",1,"a",0,"f",1);
specifyCoefficients(model2,"m",0,"d",0,"c",1,"a",0,"f",1);

applyBoundaryCondition(model1,"dirichlet","Edge",[1,2,3,4,5],"u",0);
applyBoundaryCondition(model2,"dirichlet","Edge",[1,2,3,4,5],"u",0);

u_temp1=solvepde(model1); % Pour imposer les condition aux limites des bords intérieurs
u_temp2=solvepde(model2);
u_temp1=u_temp1.NodalSolution;
u_temp2=u_temp2.NodalSolution;

cd=findNodes(model1.Mesh,"region","Edge",3);
cg=findNodes(model2.Mesh,"region","Edge",5);

u_temp1(cd)=u_sol_circle1(cd);
u_temp2(cg)=u_sol_circle2(cg);

[p1,e1,t1]=meshToPet(model1.Mesh);
F1 = pdeInterpolant(p1,t1,u_temp1);
cl_cd_1=@(region,state) evaluate(F1,region.x,region.y);

[p2,e2,t2]=meshToPet(model2.Mesh);
F2 = pdeInterpolant(p2,t2,u_temp2);
cl_cg_2=@(region,state) evaluate(F2,region.x,region.y);

applyBoundaryCondition(model1,"dirichlet","Edge",[1,2,4,5],"u",0);
applyBoundaryCondition(model1,"dirichlet","Edge",3,"u",cl_cd_1);
applyBoundaryCondition(model2,"dirichlet","Edge",[1,2,3,4],"u",0);
applyBoundaryCondition(model2,"dirichlet","Edge",5,"u",cl_cg_2);

bordmodel1=findNodes(model1.Mesh,"region","Edge",[1,2,3,4,5]);
bordmodel2=findNodes(model2.Mesh,"region","Edge",[1,2,3,4,5]);

freenodes1=setdiff(1:size(model1.Mesh.Nodes(1,:),2),bordmodel1); % Récupère les indices des noeuds libres
freenodes2=setdiff(1:size(model2.Mesh.Nodes(1,:),2),bordmodel2);

AssembleMatrix1=assembleFEMatrices(model1,"nullspace");
AssembleMatrix2=assembleFEMatrices(model2,"nullspace");

fmodel1_free=AssembleMatrix1.Kc*(AssembleMatrix1.B\(u_sol_circle1-AssembleMatrix1.ud));
fmodel2_free=AssembleMatrix2.Kc*(AssembleMatrix2.B\(u_sol_circle2-AssembleMatrix2.ud));

fmodel1=zeros(size(model1.Mesh.Nodes(1,:),2),1);
fmodel2=zeros(size(model2.Mesh.Nodes(1,:),2),1);

fmodel1(freenodes1)=fmodel1_free;
fmodel2(freenodes2)=fmodel2_free;
fmodel1(bordmodel1)=u_sol_circle1(bordmodel1); 
fmodel2(bordmodel2)=u_sol_circle2(bordmodel2);

fsource_circle1=@(location,state) function_source_f(model1,fmodel1,location,state);
fsource_circle2=@(location,state) function_source_f(model2,fmodel2,location,state);

specifyCoefficients(model1,"m",0,"d",0,"c",1,"a",0,"f",fsource_circle1);
specifyCoefficients(model2,"m",0,"d",0,"c",1,"a",0,"f",fsource_circle2);

cd=findNodes(model1.Mesh,"region","Edge",3);
cg=findNodes(model2.Mesh,"region","Edge",5);

y_bord={u_sol_circle1(cd); u_sol_circle2(cg)};

y0=zeros(size(cd));

nb_iter_schwarz=10;
eps_arret_schwarz=1e-12;
nb_cycle_aitken=5;

[res_bord1, res_mod1, err_aitkenSVD] = SchwarzAitkenSVD_2c(model1, model2, y0, y_bord, nb_iter_schwarz, eps_arret_schwarz, nb_cycle_aitken);
[res_bord2, res_mod2, err_aitken] = SchwarzAitken_2c(model1, model2, y0, y_bord, nb_iter_schwarz, eps_arret_schwarz, nb_cycle_aitken);
[cell_all_iter, cell_all_iter_bord, res_mod_gauche, res_mod_droit, err_schwarz] = iter_solve_2c(model1, model2, max(length(err_aitkenSVD),length(err_aitken)), y0, y_bord, eps_arret_schwarz);
res_mod_gauche_nodes=res_mod1{1}.NodalSolution;
res_mod_gauche_nodes(cd)=res_bord1{1};

res_mod_droit_nodes=res_mod1{2}.NodalSolution;
res_mod_droit_nodes(cg)=res_bord1{2};

% [p,e,t]=meshToPet(model1.Mesh);
% F1 = pdeInterpolant(p,t,res_bord_gauche);
% cl_cg_2=@(region,state) evaluate(F1,region.x,region.y); %condition initiale

% applyBoundaryCondition(model2,"dirichlet","Edge",[1,2,3,4],"u",0);
% applyBoundaryCondition(model2,"dirichlet","Edge",5,"u",cl_cg_2);

% results2=solvepde(model2);

figure(3)
subplot(1,2,1)
pdeplot(model1.Mesh,"XYData",res_mod_gauche_nodes) 
title("Aitken SVD (Partie Gauche)")
axis equal
subplot(1,2,2)
pdeplot(model2.Mesh,"XYData",res_mod_droit_nodes)
title("Aitken SVD (Partie Droite)")
axis equal
saveas(gcf,"Res_After_AccSVD.jpg")

figure(4)
semilogy(1:length(err_schwarz),err_schwarz,1:length(err_aitken),err_aitken,1:length(err_aitkenSVD),err_aitkenSVD)
title("Comparaison de la convergence entre Schwarz, Aitken et Aitken SVD")
legend("Schwarz","Aitken","Aitken SVD")
xlabel("Itération")
ylabel("Résidu")
saveas(gcf,"Residu.jpg")