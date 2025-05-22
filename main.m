clear all
close

R=2;
L=3;
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

figure(1)
pdegplot(model1,"EdgeLabels","on","FaceLabels","on");
title("Décomposition Face/Arête de la géométrie deux cercles")
axis equal
saveas(gcf,"geom.jpg")

cd=findNodes(model1.Mesh,"region","Edge",3);

specifyCoefficients(model1,"m",0,"d",0,"c",1,"a",0,"f",1);
specifyCoefficients(model2,"m",0,"d",0,"c",1,"a",0,"f",1);

cl_cd=0;
applyBoundaryCondition(model1,"dirichlet","Edge",[1,2,4,5],"u",0);
applyBoundaryCondition(model1,"dirichlet","Edge",3,"u",cl_cd);

results=solvepde(model1);
u=results.NodalSolution;
u_cd=u(cd);

figure(2)
pdeplot(model1.Mesh.Nodes(1,cd),model1.Mesh.Nodes(2,cd),'XYData',u_cd,'Mesh','on');
saveas(gcf,"Results.jpg")