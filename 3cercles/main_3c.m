clear 
close all

addpath("../commun")

R=2;
L=2; %Attention peut changer la numérotation des arêtes
y=0; 
x=0;

long_EF_max=0.1;
ord_EF="Quadratic";

model1=create1circle3intermodel(x,y,L,R,long_EF_max,ord_EF)
model2=create1circle3intermodel(x,y,-L,R,long_EF_max,ord_EF)
model3=create1circle3interbasmodel(x,y,L,R,long_EF_max,ord_EF)

figure(1)
hold on
pdegplot(model1,"EdgeLabels","on","FaceLabels","on");
pdegplot(model2,"EdgeLabels","on","FaceLabels","on");
pdegplot(model3,"EdgeLabels","on","FaceLabels","on");
hold off
title("Geometrie")
axis equal
saveas(gcf,"geom.jpg")

specifyCoefficients(model1,"m",0,"d",0,"c",1,"a",0,"f",1);
specifyCoefficients(model2,"m",0,"d",0,"c",1,"a",0,"f",1);
specifyCoefficients(model3,"m",0,"d",0,"c",1,"a",0,"f",1);

c1=findNodes(model1.Mesh,"region","Edge",[2,3,4]);
c2=findNodes(model2.Mesh,"region","Edge",[6,1,2]);
c3=findNodes(model3.Mesh,"region","Edge",[3,4,5]);

y0={zeros(size(c1))',zeros(size(c2))',zeros(size(c3))'};

nb_iter_schwarz=15;
eps_arret_schwarz=1e-12;
nb_cycle_aitken=10;
disp("Schwarz - Aitken SVD");
[res_bord, res_mod, err_aitkenSVD] = SchwarzAitkenSVD_3c(model1, model2, model3, y0, nb_iter_schwarz, eps_arret_schwarz, nb_cycle_aitken);
% disp("Schwarz - Aitken");
% [res_bord2, res_mod2, err_aitken] = SchwarzAitken_3c(model1, model2, model3, y0, nb_iter_schwarz, eps_arret_schwarz, nb_cycle_aitken);
disp("Schwarz, nb ité = " + length(err_aitkenSVD)+1);
[cell_all_iter, cell_all_iter_bord, res_mod_c1, res_mod_c2, res_mod_c3, err_schwarz] = iter_solve_3c(model1, model2, model3, length(err_aitkenSVD)+1, y0, eps_arret_schwarz);

res_mod_c1_nodes=res_mod{1}.NodalSolution; 
res_mod_c1_nodes(c1)=res_bord{1};

res_mod_c2_nodes=res_mod{2}.NodalSolution;
res_mod_c2_nodes(c2)=res_bord{2};

res_mod_c3_nodes=res_mod{3}.NodalSolution;
res_mod_c3_nodes(c3)=res_bord{3};

figure(2)
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

% figure(3)
% semilogy(1:length(err_schwarz),err_schwarz,1:length(err_aitken),err_aitken,1:length(err_aitkenSVD),err_aitkenSVD)
% title("Comparaison de la convergence entre Schwarz, Aitken et Aitken SVD")
% legend("Schwarz","Aitken","Aitken SVD")
% xlabel("Itération")
% ylabel("Résidu")
% saveas(gcf,"Residu.jpg")

figure(3)
semilogy(1:length(err_schwarz),err_schwarz,1:length(err_aitkenSVD),err_aitkenSVD)
title("Comparaison de la convergence entre Schwarz et Aitken SVD")
legend("Schwarz","Aitken SVD")
xlabel("Itération")
ylabel("Résidu")
saveas(gcf,"Residu.jpg")