function res=function_source_f(model,vector_f,location,state)
    [p,e,t]=meshToPet(model.Mesh);
    F2=pdeInterpolant(p,t,vector_f);
    res = evaluate(F2, location.x, location.y);
    res = res'; % PARCE QUE MATLAB NE SAIT PAS QU'ON PRIORISE TOUJOURS LA GESTION DES VECTEURS COLONNES SCREUGNEUGNEU
end