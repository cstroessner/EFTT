clear
for m = 2:4
    d = m^2; % number of discs (parameters)
    rho = 1/(4*m+2); % radius of discs
    centers = zeros(d,2);
    for s = 1:m
        for t = 1:m
            centers(s+m*(t-1),:) = [4*s-1,4*t-1];
        end
    end
    centers = rho.*centers;
    
    model = createpde();
    % Coordinates
    lowerLeft  = [0,0];
    lowerRight = [1,0];
    upperRight = [1,1];
    upperLeft =  [0,1];
    % Geometry matrix
    S = [3,4, lowerLeft(1), lowerRight(1), upperRight(1), upperLeft(1), ...
        lowerLeft(2), lowerRight(2), upperRight(2), upperLeft(2)];
    gdm = S';
    % Names
    ns = 'S';
    % Set formula
    sf = 'S';
    % Invoke decsg
    g = decsg(gdm,ns,sf');
    % Import g into model using geometryFromEdges.
    geometryFromEdges(model,g);
    generateMesh(model,'Hmax',0.1);
    applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);
    fcookie = @(x) cookiepde(x,model,centers,rho);
    switch m
        case 2
            fcookie2 = @(x) fcookie(0.75*x+1.25);
        case 3
            fcookie3 = @(x) fcookie(0.75*x+1.25);
        case 4
            fcookie4 = @(x) fcookie(0.75*x+1.25);
    end
end
save('testfunctions/cookiepde','fcookie2','fcookie3','fcookie4')

function val = cookiepde(p,model,centers,rho)
f = @(location,state) diffusion(location,state,centers,rho,p);
specifyCoefficients(model,'m',0,...
    'd',0,...
    'c',f,...
    'a',0,...
    'f',1);
results = solvepde(model);
val = mean(results.NodalSolution);
end

function a = diffusion(location,state,centers,rho,p)
x = [location.x',location.y'];
pmu = sum(p.*(pdist2(x,centers)<=rho),2);
a = ((pmu==0) .* 1 + (pmu~=0).*pmu)';
end