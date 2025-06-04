function model=create1circle3interbasmodel(x,y,L,R,s_max,order) %La définition des C1,C2,C3 est différente de create1circle3intermodel
    % x,y : Coordonnées du centre de la construction globale
    % L : écart entre les deux centres des cercles
    % R : rayons des cercles
    % s_max : longueur max des éléments finis 
    % order : Ordre de l'élément fini utilisé 

    C2=[1,x-L/2,y,R]'; % Le premier 1 précise la forme "cercle"
    C3=[1,x+L/2,y,R]';
    C1=[1,x+L/2-L/2,y-abs(L),R]';

    ens_shape=[C1,C2,C3];
    ns=char('C1', 'C2', 'C3')';
    sf="(C1-C2-C3)+(C1&C2)+(C1&C3)";
    geom=decsg(ens_shape,sf,ns); % décomposition en faces et arêtes. (EF, Table de Connectivité)

    model=createpde; 
    geometryFromEdges(model,geom); %création d'une géométrie compréhensible par PDE Toolbox

    generateMesh(model,"hmax",s_max,"GeometricOrder",order) %création du maillage
end