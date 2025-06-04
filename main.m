clear 
close all

R=2;
L=3.9;
y=0; 
x=0;

long_EF_max=0.2;
ord_EF="Quadratic";

% [model2c,mesh2c]=create2circlemesh(x,y,L,R,long_EF_max,ord_EF)

% %Nf2 = findNodes(model2c.Mesh,"region","Face",[2,3]);

% % figure(1)
% % pdegplot(model2c,"EdgeLabels","on","FaceLabels","on");
% % title("Décomposition Face/Arête de la géométrie deux cercles")
% % axis equal
% % saveas(gcf,"Geom2C.jpg")

% % figure(2)
% % pdemesh(model2c)
% % title("Maillage pour la géométrie deux cercles")
% % axis equal

% % figure(3)
% % hold on
% % pdemesh(model2c)
% % plot(model2c.Mesh.Nodes(1,Nf2),model2c.Mesh.Nodes(2,Nf2),"ok","MarkerFaceColor","g")
% % title("Appel Maillage Face 2")
% % axis equal
% % saveas(gcf, "CallFace.jpg")
% % hold off

% Er = findElements(model2c.Mesh,"region","Face",[1,3]);
% figure(4)
% pdemesh(model2c.Mesh.Nodes,model2c.Mesh.Elements(:,Er),"EdgeColor","green")
% saveas(gcf,"Mesh.jpg")

% % model1c=create1circlemodel(x-L,y,R,long_EF_max,ord_EF)
% % figure(5)
% % pdemesh(model1c)
% % saveas(gcf,"1cercle.jpg")

% specifyCoefficients(model2c,"m",0,"d",0,"c",1,"a",0,"f",0,"Face",2);
% specifyCoefficients(model2c,"m",0,"d",0,"c",1,"a",0,"f",0,"Face",3);
% specifyCoefficients(model2c,"m",0,"d",0,"c",1,"a",0,"f",1,"Face",1);
% applyBoundaryCondition(model2c,"dirichlet","Edge",1:model2c.Geometry.NumEdges,"u",0);

% results = solvepde(model2c);
% u = results.NodalSolution;

% figure(6)
% hold on 
% pdemesh(model2c)
% %pdeplot(model2c.Mesh.Nodes,model2c.Mesh.Elements(:,Er),"XYData",u)
% pdeplot(model2c,"XYData",u)
% title("Numerical Solution");
% xlabel("x")
% ylabel("y")
% hold off
% saveas(gcf,"Results.jpg")


% model1cinter=create1circleintermodel(x-L,y,R,long_EF_max,ord_EF)
% figure(7)
% pdemesh(model1cinter)
% saveas(gcf,"geom.jpg")

model1=create1circleintermodel(x,y,L,R,long_EF_max,ord_EF)
model2=create1circleintermodel(x,y,-L,R,long_EF_max,ord_EF)

specifyCoefficients(model1,"m",0,"d",0,"c",1,"a",0,"f",1);
specifyCoefficients(model2,"m",0,"d",0,"c",1,"a",0,"f",1);

% figure(1)
% pdegplot(model1,"EdgeLabels","on","FaceLabels","on");
% title("Décomposition Face/Arête de la géométrie deux cercles")
% axis equal
% saveas(gcf,"geom.jpg")

cd=findNodes(model1.Mesh,"region","Edge",3);
cg=findNodes(model2.Mesh,"region","Edge",5);

% specifyCoefficients(model1,"m",0,"d",0,"c",1,"a",0,"f",1);
% specifyCoefficients(model2,"m",0,"d",0,"c",1,"a",0,"f",1);

% cl_cd=0;
% applyBoundaryCondition(model1,"dirichlet","Edge",[1,2,4,5],"u",0);
% applyBoundaryCondition(model1,"dirichlet","Edge",3,"u",cl_cd);

% results=solvepde(model1);
% u_cd=interpolateSolution(results,model1.Mesh.Nodes(1,cd),model1.Mesh.Nodes(2,cd));

% Rest_Fun_Diric=@(region,state) evalue(u_cd, region.x, region.y)

y0=zeros(size(cd));
% [mat_all_iter,cell_all_iter_bord,res_mod_gauche,res_mod_droit]=iter_solve(model1,model2,10,y0); 
% res_bord_droit=acc_aitkenSVD(cell_all_iter_bord{1}); 
% res_bord_gauche=acc_aitkenSVD(cell_all_iter_bord{2});

% [mat_all_iter,cell_all_iter_bord,res_mod_gauche,res_mod_droit]=iter_solve(model1,model2,10,res_bord_droit); 
% res_bord_droit=acc_aitkenSVD(cell_all_iter_bord{1}); 
% res_bord_gauche=acc_aitkenSVD(cell_all_iter_bord{2});

[res_bord, res_mod, err_aitkenSVD] = SchwarzAitkenSVD(model1, model2, y0, 10, 1e-14, 30);
[res_bord2, res_mod2, err_aitken] = SchwarzAitken(model1, model2, y0, 10, 1e-14, 30);
[cell_all_iter, cell_all_iter_bord, res_mod_gauche, res_mod_droit, err_schwarz] = iter_solve(model1, model2, length(err_aitkenSVD)+1, y0);

% disp([res_bord_droit,cell_all_iter_bord{1}(:,end)])
% disp([res_bord_gauche,cell_all_iter_bord{2}(:,end)])
% disp(norm(res_bord_droit-cell_all_iter_bord{1}(:,end)))
% disp(norm(res_bord_gauche-cell_all_iter_bord{2}(:,end)))

res_mod_gauche_nodes=res_mod{1}.NodalSolution; 
res_mod_gauche_nodes(cd)=res_bord{1};

% [p,e,t]=meshToPet(model2.Mesh);
% F2 = pdeInterpolant(p,t,res_mod_gauche_nodes);
% cl_cd_1=@(region,state) evaluate(F2,region.x,region.y); %condition initiale

% applyBoundaryCondition(model1,"dirichlet","Edge",[1,2,4,5],"u",0);
% applyBoundaryCondition(model1,"dirichlet","Edge",3,"u",cl_cd_1);

% results1=solvepde(model1);

res_mod_droit_nodes=res_mod{2}.NodalSolution;
res_mod_droit_nodes(cg)=res_bord{2};

% [p,e,t]=meshToPet(model1.Mesh);
% F1 = pdeInterpolant(p,t,res_bord_gauche);
% cl_cg_2=@(region,state) evaluate(F1,region.x,region.y); %condition initiale

% applyBoundaryCondition(model2,"dirichlet","Edge",[1,2,3,4],"u",0);
% applyBoundaryCondition(model2,"dirichlet","Edge",5,"u",cl_cg_2);

% results2=solvepde(model2);

disp(norm(res_bord{1}-res_bord2{1}))
disp(norm(res_bord{2}-res_bord2{2}))
disp(norm(res_bord{1}-res_mod_gauche.NodalSolution(cd)))
disp(norm(res_bord{2}-res_mod_droit.NodalSolution(cg)))
disp(norm(res_bord2{1}-res_mod_gauche.NodalSolution(cd)))
disp(norm(res_bord2{2}-res_mod_droit.NodalSolution(cg)))

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