function fun=circle2_solution_manufactured(R,L)
    % fun1=@(x,y) ((x-L/2).^2 + y.^2-R^2<0).*(sin(((x-L/2).^2 + y.^2)*pi/R^2));
    % fun2=@(x,y) ((x+L/2).^2 + y.^2-R^2<0).*(sin(((x+L/2).^2 + y.^2)*pi/R^2));
    fun1=@(x,y) (1-((x-L/2).^2+y.^2)/R^2).^3.*((x-L/2).^2 + y.^2<R^2); 
    fun2=@(x,y) (1-((x+L/2).^2+y.^2)/R^2).^3.*((x+L/2).^2 + y.^2<R^2);
    fun=@(x,y) (fun1(x,y)+fun2(x,y))';
end