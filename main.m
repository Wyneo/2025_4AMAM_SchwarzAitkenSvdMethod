R=2;
L=3;
y=0; 
x=0;

long_EF_max=0.2;
ord_EF="Quadratic";

[model2c,mesh2c]=create2circlemesh(x,y,L,R,long_EF_max,ord_EF)

Nf2 = findNodes(mesh2c,"region","Face",3);


figure(1)
pdegplot(model2c,"EdgeLabels","on","FaceLabels","on");
title("Décomposition Face/Arête de la géométrie deux cercles")
axis equal

saveas(gcf,"Geom2C.jpg")

figure(2)
pdemesh(model2c)
title("Maillage pour la géométrie deux cercles")
axis equal

figure(3)
plot(mesh2c.Nodes(1,Nf2),mesh2c.Nodes(2,Nf2),"ok","MarkerFaceColor","g")
title("Appel Maillage Face 2")
axis equal
saveas(gcf, "CallFace.jpg")




