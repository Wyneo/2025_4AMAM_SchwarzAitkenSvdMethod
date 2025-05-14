R=2;
L=3;
y=0; 
x=0;

long_EF_max=0.2;
ord_EF="Quadratic";

model=create2circlemodel(x,y,L,R,long_EF_max,ord_EF)

[p,e,t] = meshToPet(model.Mesh); %decomposition PET
% p : Matrice des points
% e : Matrice des arêtes 
% t : Matrice des triangles 

nodes = p'; %extraction des noeuds
nn = size(nodes,1); %nombre de noeus

% Identification des sous-domaines (A changer vu comment c'est pas fou)
inA = vecnorm(nodes - [x+L/2,y], 2, 2) < R; %identifie les différences aux centres et en tire un booléen
inB = vecnorm(nodes - [x-L/2,y], 2, 2) < R;

idF1 = find(inA & ~inB);
idF2 = find(inA &  inB);
idF3 = find(inB & ~inA);

idF12 = union(idF1, idF2); %union des sous domaines
idF23 = union(idF2, idF3);

applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);

specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',1);

% Assemblage global
FEM = assembleFEMatrices(model); %all odes matrix
K = FEM.K; %stiffness matrix (matrice du laplacien)
F = FEM.F; %vecteur b intégrale

u = zeros(nn,1); %??
Niter = 20;

for k = 1:Niter
    I12 = idF12;
    u12 = zeros(nn,1);
    K12 = K(I12,I12);
    F12 = F(I12);
    u12(I12) = K12 \ F12;

    u(I12) = u12(I12); % mise à jour partielle

    I23 = idF23;
    u23 = zeros(nn,1);
    K23 = K(I23,I23);
    F23 = F(I23);
    u23(I23) = K23 \ F23;

    u(I23) = u23(I23);
end

pdeplot(model,'XYData',u,'Mesh','on');
title('Solution FEM avec DDM sur deux cercles intersectés');
axis equal
saveas(gcf,"test.jpg")