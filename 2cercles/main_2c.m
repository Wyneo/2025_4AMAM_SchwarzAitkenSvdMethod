clear 
close all

addpath("../commun")

R=2;
L=2.8; %Attention si L>2.8 la numérotation des arêtes change
y=0; 
x=0;

long_EF_max=0.2;
ord_EF="Quadratic";

model1=create1circleintermodel(x,y,L,R,long_EF_max,ord_EF);
model2=create1circleintermodel(x,y,-L,R,long_EF_max,ord_EF);

figure(1)
hold on
pdegplot(model1,"EdgeLabels","on","FaceLabels","on");
pdegplot(model2,"EdgeLabels","on","FaceLabels","on");
hold off
title("Geometrie")
axis equal
saveas(gcf,"geom.jpg")

specifyCoefficients(model1,"m",0,"d",0,"c",1,"a",0,"f",1);
specifyCoefficients(model2,"m",0,"d",0,"c",1,"a",0,"f",1);

cd=findNodes(model1.Mesh,"region","Edge",[2,3]);
cg=findNodes(model2.Mesh,"region","Edge",[1,5]);

y0=zeros(size(cd));

nb_iter_schwarz=10;
eps_arret_schwarz=1e-12;
nb_cycle_aitken=10;
disp("Schwarz - Aitken SVD")
[res_bord1, res_mod1, err_aitkenSVD] = SchwarzAitkenSVD_2c(model1, model2, y0, nb_iter_schwarz, eps_arret_schwarz, nb_cycle_aitken);
disp("Schwarz - Aitken")
[res_bord2, res_mod2, err_aitken] = SchwarzAitken_2c(model1, model2, y0, nb_iter_schwarz, eps_arret_schwarz, nb_cycle_aitken);
disp("Schwarz, nb ité = "+length(err_aitkenSVD)+1)
[cell_all_iter, cell_all_iter_bord, res_mod_gauche, res_mod_droit, err_schwarz] = iter_solve_2c(model1, model2, length(err_aitkenSVD)+1, y0, eps_arret_schwarz);

res_mod_gauche_nodes=res_mod1{1}.NodalSolution;
res_mod_gauche_nodes(cd)=res_bord1{1};
res_mod_droit_nodes=res_mod1{2}.NodalSolution;
res_mod_droit_nodes(cg)=res_bord1{2};

figure(2)
subplot(1,2,1)
pdeplot(model1.Mesh,"XYData",res_mod_gauche_nodes) 
title("Aitken SVD (Partie Gauche)")
axis equal
subplot(1,2,2)
pdeplot(model2.Mesh,"XYData",res_mod_droit_nodes)
title("Aitken SVD (Partie Droite)")
axis equal
saveas(gcf,"Res_After_AccSVD.jpg")

figure(3)
semilogy(1:length(err_schwarz),err_schwarz,1:length(err_aitken),err_aitken,1:length(err_aitkenSVD),err_aitkenSVD)
title("Comparaison de la convergence entre Schwarz, Aitken et Aitken SVD")
legend("Schwarz","Aitken","Aitken SVD")
xlabel("Itération")
ylabel("Résidu")
saveas(gcf,"Residu.jpg")