function [d,D] = searchdirection(g,H,i,tr,algorithm)
%Pass the algorithm as a string and compare it with the name of the
%optsolver argument.
%g is the gradient of the function at the current iterate.
%Here I've passed the Hessian as per SR1 and 
if contains(algorithm,'steepest')
    d = -g; %for the steepestdescent direction is the negative of the gradient.
    D = d'*g;
elseif contains(algorithm,'newton')
    %In newton method we need hessian as well as gradient for finding the
    %direction but our hessian needs to  be postitve definite because its
    %in the denominator
    xi = 1e-4;
    n = length(g);
    while min(eig(H))<0
        H = H + xi*eye(n);
        xi = 10*xi;
        
    end
    % direction in the newton
    d = H\(-g);
    D = d'*g;
% CG algorithm is used to solve the subproblem of the Trust Region
% approximation here we solve the AX = b equation for X, x is the descent direction of the function. 
elseif contains(algorithm,'cg')
    A = H; % Initialize A as hessian 
    b= -g; % Initialize b as negative of the gradient
    X = zeros(size(g)); % X is our descent direction and we intialize as zero in CG 
    r= A'*X-b; %residual 
    p=-r;
    k= 1;
    while k<i.cgmaxiter
        if p'*A*p<0
            c = X'*X - tr^2;
            b = 2*p'*X;
            a = p'*p;
            alpha = max(roots([a b c])); %Solve the quadratic equation to find the value of the alpha
            X = X+alpha*p; % return the value of X.
            break
        else 
            alpha = (r'*r)/(p'*A*p);
        end
        if (norm(X+alpha*p)) > tr
            c = X'*X - tr^2;
            b = 2*p'*X;
            a = p'*p;
            alpha = max(roots([a b c])); %Solve the quadratic equation to find the value of the alpha
            X = X + alpha*p; % return the value of X.
            break
        else
            X = X + alpha*p; % retunr the value of X.
            r_new = r + alpha*A*p;
        end
        if norm(r_new) < i.cgopttol
            break
            
        else 
            beta = (r_new'*r_new)/(r'*r);
            p = -r_new + beta*p;
            r = r_new;
        end
        k = k+1;
    end
    d = X;  %X is the descent direction 
    D = d'*g;
% In the BFGS instead of using the actual hessian we use hessian that is
% calculated by BFGS algorithm the direction is same as that of newton only
% hessian calculation changes 
elseif contains(algorithm,'bfgs')
    xi = 1e-4;
    n = length(g);
    while min(eig(H))<0
        H = H + xi*eye(n);
        xi = 10*xi;
        
    end
    % direction in the newton
    d = H\(-g);
    D = d'*g;
end


   

     

    





    









