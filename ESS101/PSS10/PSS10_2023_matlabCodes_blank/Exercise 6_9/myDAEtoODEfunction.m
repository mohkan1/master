function [dxdt,g,z,detdg] = myDAEtoODEfunction( t, x, z,uoft, tol, verbose )
    %% This function evaluates dx/dt for a given state x and input u
    
    if verbose
        
        disp(repmat('-',[1,100]))
        fprintf('integration time t = %2.14f\n',t)
        fprintf('Running newton method\n\nit\t|g|\t\t|dz|\t\talpha\t\tdet(dg)\n')
    end
    
    niter    = 0;
    g        = inf;
    maxIter  = 30;
    lsStatus = 1;
    
    while niter <= maxIter && norm(g)> tol && lsStatus == 1
        niter = niter + 1;
        [g,dg] = algebraicRhs(x,z,uoft);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dz     = ... % your code here!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [alpha, g,lsStatus] = LineSearch(g,x,z,uoft,dz,tol);

        z = z+alpha*dz;  
        if verbose
           fprintf('%2i\t%1.1e\t\t%1.1e\t\t%1.1e\t\t%1.1e\n',niter,norm(g),norm(dz),alpha,det(dg)) 
        end
    end

    
    dxdt = differentialRhs(x,z,uoft);
    detdg = det(dg);
   
    if niter >= maxIter 
        if lsStatus == -1
            error('Linesearch failed. Details: \n\nalpha\t\t = %1.1e\n|g|\t\t = %1.1e\ndet(jacobian(g)) = %2.2e\n\n\n\n',alpha,norm(g),det(dg))
        else
            error('Max iterations in newton reached without achieveing set tolerance. Details: \n\nalpha\t\t = %1.1e\n|g|\t\t = %1.1e\ndet(jacobian(g)) = %2.2e\n\n\n\n',alpha,norm(g),det(dg))
        end
    elseif lsStatus == -1
        warning('Linesearch failed. Details: \n\nalpha\t\t = %1.1e\n|g|\t\t = %1.1e\ndet(jacobian(g)) = %2.2e\n\n\n\n',alpha,norm(g),det(dg))    
        if isnan(dz)
            error('!')
        end
    elseif isnan(dz)
        error('Direction dz is NaN. Details: \n\nalpha\t\t = %1.1e\n|g|\t\t = %1.1e\ndet(jacobian(g)) = %2.2e\n\n\n\n',alpha,norm(g),det(dg))    
    else
        if verbose
        fprintf('\nSolution found\n')
        end
    end
    
end

function [alpha, gnew,lsStatus] = LineSearch(g,x,z,uoft,dz,tolerance)

alpha    = 1;
LS       = true;
lsStatus = 1;
while LS
    gnew = algebraicRhs(x,z+alpha*dz,uoft);
    
    if (norm(gnew) < norm(g)) || (norm(gnew)<tolerance) % accept step
        LS = false;
    else
        alpha = 0.7*alpha;
        if alpha < 1e-14
            lsStatus = -1;
            break
        end
    end
end
end

