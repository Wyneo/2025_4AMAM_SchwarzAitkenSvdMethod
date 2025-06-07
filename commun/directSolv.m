function res=directSolv(model,title_name,file_name)
    res = solvepde(model);
    res_plot = res.NodalSolution; % solution at the point of the mesh
    figure
    pdeplot(model,"XYdata",res_plot)
    title(title_name)
    axis equal
    saveas(gcf,file_name)
end