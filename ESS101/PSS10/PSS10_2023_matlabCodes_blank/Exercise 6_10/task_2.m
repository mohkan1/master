%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%               PSS on Solution to DAE, fully implicit DAE Case
%                   
% Written by Robert Hult, adapted from Sebastien Gros.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear


% WRITE FUNCTIONS ---------------------------------------------------------
subAssignment = 'a';
switch subAssignment
    case 'a'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % z  = ... % your code here
        % x  = ... % your code here
        % dx = ... % your code here
        % u  = ... % your code here
        % 
        % F = ... % your code here
        %%%%%%%%%%%%%%%%%%%%%%%%%%
    
        w      = [dx;z]; 
        dF     = jacobian(F,w); 
        x0     = 2;
        tstop  = 10; % termination time of simulation

        un = @(t) 1;
    case 'b' % helicopter from assignment 1
           % try this out!
    otherwise
        warning('non-existing subassignment')
        return
end
matlabFunction(F,dF,'File','implicitRHS','vars',{w,x,u});




% SIMULATE ODE ------------------------------------------------------------
tolerance = 1e-6;
printNewtonProgress = true;
z0  = zeros(length(z),1);   % initial guess for newton
dx0 = zeros(length(dx),1); % initial guess for newton
odeWrapper = @(t,x) myImplicitDAEtoODEfunction( t,dx0,x,z0,un(t),tolerance, printNewtonProgress );

tstart  = 0;
odesol  = ode45(odeWrapper,[tstart,tstop],x0);



% PLOT RESULTS ------------------------------------------------------------
% PLOT X
figure(1)
clf
nx= length(x);

for i = 1:nx
   subplot(nx+1,1,i)
   hold on
   grid on
   plot(odesol.x,odesol.y(i,:))
   xlabel('$t$','interpreter','latex')
   ylabel(sprintf('$x_%i$',i),'interpreter','latex')
end
% plot |F|
for i = 1:length(odesol.x)
    t = odesol.x(i);
    x = odesol.y(:,i);
    [~, Ftmp]= myImplicitDAEtoODEfunction( t,dx0,x,z0,un(t),tolerance, false );
    Fnorm(i) = norm(Ftmp);
end
subplot(nx+1,1,nx+1)
semilogy(odesol.x,Fnorm)
hold on
grid on

xlabel('t')
ylabel('$|F(\dot x,x,z,u)|$','interpreter','latex')



