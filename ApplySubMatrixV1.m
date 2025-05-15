R=2;
L=3;
y=0; 
x=0;

long_EF_max=0.2;
ord_EF="Quadratic";

model=create2circlemodel(x,y,L,R,long_EF_max,ord_EF)

figure(1)
pdegplot(model,"EdgeLabels","on","FaceLabels","on");
title("Décomposition Face/Arête de la géométrie deux cercles")
axis equal
saveas(gcf,"Geom2C.jpg")

nn = size(model.Mesh.Nodes,2); %nombre de noeuds

% Identification des sous-domaines 
idF23=transpose(findNodes(model.Mesh,"region","Face",[2,3]));
idF23_bord=findNodes(model.Mesh,"region","Edge",[1:5]);
idF12=transpose(findNodes(model.Mesh,"region","Face",[1,3]));
idF12_bord=findNodes(model.Mesh,"region","Edge",[5:10]);
idF23=idF23(~ismember(idF23, idF23_bord));
idF12=idF12(~ismember(idF12, idF12_bord));

applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);

specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',1);

% Assemblage global
FEM = assembleFEMatrices(model); %all odes matrix
K = FEM.K; %stiffness matrix (matrice du laplacien)
F = FEM.F; %vecteur b intégrale

u = zeros(nn,1);
Niter = 20;

for k = 1:Niter
    I12 = idF12;
    u12 = zeros(nn,1);
    K12 = K(I12,I12);
    F12 = F(I12);
    u12(I12) = K12 \ F12;

    u(I12) = u12(I12); % mise à jour partielle

    %applyBoundaryCondition(model,'dirichlet','Edge',10,'u',u(arc_cercle_gauche),"InternalBC",true);

    I23 = idF23;
    u23 = zeros(nn,1);
    K23 = K(I23,I23);
    F23 = F(I23);
    u23(I23) = K23 \ F23;

    u(I23) = u23(I23);
    %applyBoundaryCondition(model,'dirichlet','Edge',3,'u',u(arc_cercle_droite),"InternalBC",true);
end

pdeplot(model,'XYData',u,'Mesh','on');
title('Solution FEM avec DDM sur deux cercles intersectés');
axis equal
saveas(gcf,"test.jpg")

