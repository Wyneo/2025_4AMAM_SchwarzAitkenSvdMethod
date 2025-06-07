clear 
close all

addpath("./commun")

R=2;
L=2.5; 
y=0; 
x=0;

circle2_sol=circle2_solution_manufactured(R,L);

long_EF_max=1;
ord_EF="Linear";

model1=create1circleintermodel(x,y,L,R,long_EF_max,ord_EF);
model2=create1circleintermodel(x,y,-L,R,long_EF_max,ord_EF);

u_fun=circle2_solution_manufactured(R,L);

% ! Tentative d'appliquer l'opérateur de Laplace Discret, puis d'utiliser la solution pour définir la source

u_sol_circle1=u_fun(model1.Mesh.Nodes(1,:),model1.Mesh.Nodes(2,:));
u_sol_circle2=u_fun(model2.Mesh.Nodes(1,:),model2.Mesh.Nodes(2,:));

specifyCoefficients(model1,"m",0,"d",0,"c",1,"a",0,"f",1);
specifyCoefficients(model2,"m",0,"d",0,"c",1,"a",0,"f",1);

applyBoundaryCondition(model1,"dirichlet","Edge",[1,2,3,4,5],"u",0);
applyBoundaryCondition(model2,"dirichlet","Edge",[1,2,3,4,5],"u",0);

u_temp1=solvepde(model1);
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

stiffmatrix1=assembleFEMatrices(model1,"K").K;
stiffmatrix2=assembleFEMatrices(model2,"K").K;

[P1,R1,C1]=equilibrate(stiffmatrix1);
[P2,R2,C2]=equilibrate(stiffmatrix2);
stiffmatrix1=R1*P1*stiffmatrix1*C1;
stiffmatrix2=R2*P2*stiffmatrix2*C2;

fmodel1=stiffmatrix1\u_sol_circle1';
fmodel2=stiffmatrix2\u_sol_circle2';

disp(norm(stiffmatrix1*fmodel1-u_sol_circle1',Inf));
disp(norm(stiffmatrix2*fmodel2-u_sol_circle2',Inf));

fsource_circle1=@(location,state) function_source_f(model1,fmodel1,location,state);
fsource_circle2=@(location,state) function_source_f(model2,fmodel2,location,state);

specifyCoefficients(model1,"m",0,"d",0,"c",1,"a",0,"f",fsource_circle1);
specifyCoefficients(model2,"m",0,"d",0,"c",1,"a",0,"f",fsource_circle2);

cd=findNodes(model1.Mesh,"region","Edge",3);
cg=findNodes(model2.Mesh,"region","Edge",5);

y0=zeros(size(cd));

[res_bord, res_mod, err_aitkenSVD] = SchwarzAitkenSVD_2c(model1, model2, y0, 10, 1e-12, 30);
[res_bord2, res_mod2, err_aitken] = SchwarzAitken_2c(model1, model2, y0, 10, 1e-12, 30);
[cell_all_iter, cell_all_iter_bord, res_mod_gauche, res_mod_droit, err_schwarz] = iter_solve2c(model1, model2, max(length(err_aitkenSVD),length(err_aitken))+1, y0, 1e-12);
res_mod_gauche_nodes=res_mod{1}.NodalSolution; 
res_mod_gauche_nodes(cd)=res_bord{1};

res_mod_droit_nodes=res_mod{2}.NodalSolution;
res_mod_droit_nodes(cg)=res_bord{2};

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