R=2;
L=3;
y=0; 
x=0;

long_EF_max=0.2;
ord_EF="Quadratic";

model2c=create2circlemesh(x,y,L,R,long_EF_max,ord_EF)

figure(1)
pdegplot(model2c,"EdgeLabels","on","FaceLabels","on");
title("Décomposition Face/Arête de la géométrie deux cercles")
axis equal

saveas(gcf,"Geom2C.jpg")

figure(2)
pdemesh(model2c)
title("Maillage pour la géométrie deux cercles")
axis equal

saveas(gcf, "Maillage2C.jpg")


