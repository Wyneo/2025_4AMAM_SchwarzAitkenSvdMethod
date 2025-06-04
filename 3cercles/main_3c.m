clear 
close all

R=2;
L=3;
y=0; 
x=0;

long_EF_max=0.05;
ord_EF="Quadratic";

model1=create1circle3intermodel(x,y,L,R,long_EF_max,ord_EF)
model2=create1circle3intermodel(x,y,-L,R,long_EF_max,ord_EF)
model3=create1circle3interbasmodel(x,y,L,R,long_EF_max,ord_EF)

specifyCoefficients(model1,"m",0,"d",0,"c",1,"a",0,"f",1);
specifyCoefficients(model2,"m",0,"d",0,"c",1,"a",0,"f",1);
specifyCoefficients(model3,"m",0,"d",0,"c",1,"a",0,"f",1);

figure(1)
hold on
pdegplot(model1,"EdgeLabels","on","FaceLabels","on");
%pdegplot(model2,"EdgeLabels","on","FaceLabels","on");
%pdegplot(model3,"EdgeLabels","on","FaceLabels","on");
hold off
title("Décomposition Face/Arête de la géométrie deux cercles")
axis equal
saveas(gcf,"geom.jpg")

c1=findNodes(model1.Mesh,"region","Edge",[2,3,4]);
c2=findNodes(model2.Mesh,"region","Edge",[6,1,2]);
c3=findNodes(model3.Mesh,"region","Edge",[3,4,5]);

y0={zeros(size(c1)),zeros(size(c2)),zeros(size(c3))};

[res_bord, res_mod, list_residu] = SchwarzAitken_3c(model1, model2, model3, y0, 10, 1e-6, 5);

res_mod_c1_nodes=res_mod{1}.NodalSolution; 
res_mod_c1_nodes(c1)=res_bord{1};

res_mod_c2_nodes=res_mod{2}.NodalSolution;
res_mod_c2_nodes(c2)=res_bord{2};

res_mod_c3_nodes=res_mod{3}.NodalSolution;
res_mod_c3_nodes(c3)=res_bord{3};

figure(3)
subplot(1,3,1)
pdeplot(model1.Mesh,"XYData",res_mod_c1_nodes) 
title("Aitken SVD (Cercle 1)")
axis equal
subplot(1,3,2)
pdeplot(model2.Mesh,"XYData",res_mod_c2_nodes)
title("Aitken SVD (Cercle 2)")
axis equal
subplot(1,3,3)
pdeplot(model3.Mesh,"XYData",res_mod_c3_nodes)
title("Aitken SVD (Cercle 3)")
axis equal
saveas(gcf,"Res_After_AccSVD.jpg")

figure(4)
subplot(3,1,1)
plot(1:length(list_residu{1}),list_residu{1})
title("Résidu de convergence (Intersection 1)")
xlabel("Itération")
ylabel("Résidu")
subplot(3,1,2)
plot(1:length(list_residu{2}),list_residu{2})
title("Résidu de convergence (Intersection 2)")
xlabel("Itération")
ylabel("Résidu")
subplot(3,1,3)
plot(1:length(list_residu{3}),list_residu{3})
title("Résidu de convergence (Intersection 3)")
xlabel("Itération")
ylabel("Résidu")
saveas(gcf,"Residu.jpg")