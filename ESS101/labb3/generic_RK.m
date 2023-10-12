function xkp1 = generic_RK(butcher,xk,dt,f,tk,u)

[s, ~] = size(butcher.a);

k = zeros(s,length(xk));
    
bsum = zeros(1,length(xk));

for i=1:s
    asum = zeros(1,length(xk));
    
    for j=1:s
        asum = asum + butcher.a(i,j).*k(j,:);
    end
    try
        k(i,:) = f(xk + dt.*asum,u(tk+butcher.c(i)*dt));
    catch 
        k(i,:) = f(xk + dt.*asum);
    end
        bsum = bsum + butcher.b(i).*k(i,:);
end
xkp1 = xk + dt.*bsum;
end
