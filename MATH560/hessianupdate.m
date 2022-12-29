function H = hessianupdate(H,p,x,a,d,g,algorithm,i)
% in the hessian update I've divided the function in three differnent
% parts, if the algorithm is SR1 then update the hessian as per SR1
% algorithm if it is BFGS then calculate as per BFGS algorithm, other than
% those two calculate as per feval(p,x,2)
s = a*d; % s = x(k+1) - x(k)
y = feval(p,x,1) - feval(p,x-a*d,1); %y = g(xk+1)-g(xk)
if a == 0  %if we don't accept the radius of the trust region then we don't update the hessian
    return
end
if contains(algorithm,'sr1')
    if a ~= 0
        if abs(s'*(y - H*s)) >= i.sr1updatetol*norm(s)*norm(y-H*s) %SR1 update tolerance condition check
            H = H + ((y - H*s)*(y-H*s)')/((y-H*s)'*s);             %update the hessian by SR1 algorithm if above condition is true
        end
    end
%Hessian update for BFGS
elseif contains(algorithm,'bfgs')
    xi = 1e-4;
    n = length(g);
    while min(eig(H))<0
        H = H + xi*eye(n); %from preventing hessian to become zero matrix
        xi = 10*xi;
        
    end
    % I followed the dampened BFGS algorithm from book N&W
    if s'*y >= 0.2*s'*H*s
        theta = 1;
    end
    if s'*y < 0.2*s'*H*s
        theta = (0.8*s'*H*s)/(s'*H*s - s'*y);
    end
    r = theta*y + (1-theta)*H*s;
    H = H - ((H*s*s'*H))/((s'*H*s)) + ((r*r')/(s'*r)); %update the hessian 

else 
    H = feval(p,x,2);   %calculate the hessian for other algorithm
end
end

