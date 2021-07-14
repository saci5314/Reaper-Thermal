function M = mach(eng) 
    %%% This function finds mach numbers along nozzle 
    
    r = eng.r;
    gamma = eng.gamma;
    A_t = pi*min(r)^2;
    
    M = zeros(1, length(r));

    gp1 = gamma+1; gm1 = gamma-1;
    a = .5*(gp1)/(gm1);

    search = [1+1e-6 30];

    for j = 1:length(r)
        A = pi*r(stat)^2;

        if r(j) == min(r)
            search = [1e-6 1];
        elseif r(j) == r(j-1)
            for i = j:length(r)
                M(j:length(r)) = linspace(M(j-1), 0, length(r)-J);
            end
            break
        end

        problem.objective = @(M) (A_t/A)*(1+gm1*M^2/2)^a-M*(gp1/2)^a;
        problem.solver = 'fzero';
        problem.options = optimset(@fzero);
        problem.x0 = search;
        M(j) = fzero(problem);
    end
end