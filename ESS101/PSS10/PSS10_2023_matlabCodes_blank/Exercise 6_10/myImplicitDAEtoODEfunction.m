function [dxdt,F] = myImplicitDAEtoODEfunction( t, dx,x, z,uoft, tol, verbose )
    %% This function evaluates dx/dt for a given state x and input u
    
    if verbose
        display(['------------------- t = ',num2str(t),' ------------------------------------'])
        
        fprintf('Running newton method\n\nit\t|F|\t\t|dw|\t\talpha\t\tdet(dF)\n')
    end
    
    niter    = 0;
    F        = inf;
    maxIter  = 30;
    lsStatus = 1;
    
    w = [dx;z]; % initial guess
    
    while niter <= maxIter && norm(F)> tol && lsStatus == 1
        niter = niter + 1;
        [F,dF] = implicitRHS(w,x,uoft);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dw     = ... % your code here
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [alpha, F,lsStatus] = LineSearch(F,x,uoft,w,dw,tol);

        w = w+alpha*dw;  
        if verbose
           fprintf('%2i\t%1.1e\t\t%1.1e\t\t%1.1e\t\t%1.1e\n',niter,norm(F),norm(dw),alpha,det(dF)) 
        end
    end

    
    dxdt = w(1:length(dx));
    
   
    if niter >= maxIter 
        error('Max iterations in newton reached without achieveing set tolerance')
    elseif lsStatus == -1
        keyboard
        error('Linesearch failed. Details: \n\nalpha\t\t = %1.1e\n|g|\t\t = %1.1e\ndet(jacobian(F)) = %2.2e\n\n\n\n',alpha,norm(F),det(dF))    
    else
        if verbose
        fprintf('\nSolution found\n')
        end
    end
end

function [alpha, Fnew,lsStatus] = LineSearch(F,x,uoft,w,dw,tolerance)

alpha    = 1;
LS       = true;
lsStatus = 1;
while LS
    Fnew = implicitRHS(w+dw*alpha,x,uoft);
    
    if (norm(Fnew) < norm(F)) || (norm(Fnew)<tolerance) % accept step
        LS = false;
    else
        alpha = 0.7*alpha;
        if alpha < tolerance/10
            lsStatus = -1;
            break
        end
    end
end
end

