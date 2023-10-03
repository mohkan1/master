%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%               PSS on Solution to DAE, semi-explict DAE Case
%                   
% Written by Robert Hult, adapted from Sebastien Gros.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

% WRITE FUNCTIONS ---------------------------------------------------------
subAssignment = 'b';
switch subAssignment
    case 'b'
        z = sym('z',[2,1],'real');
        x = sym('x',[2,1],'real');
        u = sym('u',[1,1],'real');

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % f = ... % your code here
        % g = ... % your code here
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Dg     = jacobian(g,z); 
        x0     = [4; 2];
        tstop  =  10; % termination time of simulation

        un = @(t) 0;
    case 'c'
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % z = ... % your code here
        % x = ... % your code here
        % u = ... % your code here
        % f = ... % your code here
        % g = ... % your code here
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Dg     = jacobian(g,z);
        x0     = [1; 0];
        tstop  =  10; % termination time of simulation
        un     = @(t) 0;
        

    case 'd'
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % z = ... % your code here
        % x = ... % your code here
        % u = ... % your code here
        % f = ... % your code here
        % g = ... % your code here
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Dg     = jacobian(g,z);
        x0     = [4; 2; 0];
        tstop  =  10; % termination time of simulation
        un = @(t) 0;
    otherwise
        warning('non-existing subassignment')
end
matlabFunction(g,Dg,'File','algebraicRhs','vars',{x,z,u});
matlabFunction(f,'File','differentialRhs','vars',{x,z,u});


% SIMULATE ODE ------------------------------------------------------------
tolerance = 1e-8;
printNewtonProgress = true;
z0 = zeros(length(z),1); % initial guess for newton
odeWrapper = @(t,x) myDAEtoODEfunction( t,x,z0,un(t),tolerance, printNewtonProgress );

tstart = 0;
odesol  = ode45(odeWrapper,[tstart,tstop],x0,odeset('RelTol',1e-8,'AbsTol',1e-8));



% PLOT RESULTS ------------------------------------------------------------
% PLOT X
figure(1)
clf
nx= length(x);

for i = 1:nx
   subplot(nx+2,1,i)
   hold on
   grid on
   plot(odesol.x,odesol.y(i,:))
   xlabel('t')
   ylabel(sprintf('x_%i',i))
end
% plot |g|
for i = 1:length(odesol.x)
    t = odesol.x(i);
    x = odesol.y(:,i);
    [~, gtmp,z,detdg(i)]= myDAEtoODEfunction( t, x, z0,un(t),tolerance, false );
    
    gnorm(i) = norm(gtmp); 
end
subplot(nx+2,1,nx+1)
semilogy(odesol.x,gnorm)
hold on
grid on
xlabel('t')
ylabel('$|g(x,z)|$','interpreter','latex')

subplot(nx+2,1,nx+2)
semilogy(odesol.x,abs(detdg))
hold on
grid on
xlabel('t')
ylabel('$\mathrm{det}\left(\frac{\partial g }{\partial z}\right)$','interpreter','latex')




