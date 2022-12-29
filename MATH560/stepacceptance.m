function [a,f,tr] = stepacceptance(p,x,f,d,D,i,H,tr,algorithm)

% this function gives the stepsize by different input of algorithms  such
% as backtrack, wofe and trustregion.

if contains(algorithm,'backtrack')  %backtracking algorithm
    c = 0.8;    %intializing the parameter
    a= 1;       %taking the initialstepsize as 1
    rho = 0.8;  %initializing the parameter
    
    while  feval(p,x+a*d,0) > feval(p,x,0)+c*a*feval(p,x,1)'*d
        a = rho*a;
    end
elseif contains(algorithm,'wolfe') %This one is wolfe algorithm from the N&W book
   %Here we apply wolfe algorithm  
   alpha_0 = 0;     %Initializing the alpha value 
   alpha_max = 5;   
   alpha_prev = alpha_0;
   alpha_i = (alpha_0+alpha_max)/2;
   k = 1;
    while k < 1e+2
   
       if feval(p,x+alpha_i*d,0) > f +i.c1ls*alpha_i*feval(p,x,1)'*d || ((feval(p,x+alpha_i*d,0) >= feval(p,x+alpha_prev*d,0)) && i > 1 )
           a = zoom(alpha_prev,alpha_i,p,x,d,i); %calling a zoom function, zoom function is defined in the last after the function definition of the stepacceptance ends
           return;                      %Zoom function gives the best stepsize for our function in the descent direction. 
       end
       phid = feval(p,x+alpha_i*d,1)'*d;
       if abs(phid) <= -i.c2ls*feval(p,x,1)'*d
           a = alpha_i;             % giving the best stepsize to use.
           return;
       end
       if phid >= 0
           a = zoom(alpha_next,alpha_prev,p,x,d,i); %zoom gives the best stepsize to use
           k = k+1;
       end
       alpha_i = (alpha_i+alpha_max)/2;
      
    end 
% Trust region algoritm is below I followed it from the slides 
elseif contains(algorithm,'trustregion')
    m_k = @(d) f + feval(p,x,1)'*d + 0.5*d'*H*d;
    m = m_k(d);
    
    rho_k = @(d) ((f - feval(p,x+d,0))/(f-m));
    rho = rho_k(d);
    gamma = 2;    
    if rho >= i.c2tr % accept the stepsize 
        tr = 2*tr; %increase the trust region size
        a = 1;
    
    elseif rho>i.c1tr && rho<i.c2tr % We accept the stepsize 
        a = 1;                      % We maintain the radius no change in the radius size
    elseif rho <= i.c1tr            % We reduce the trust region and do not accept the stepsize
        tr = tr/gamma;
        a= 0;


    end
end

end
% Zoom function definition is here I followed it from the N&W book
function alpha = zoom(alpha_prev,alpha_next,p,x,d,i)
        alpha_hi = alpha_next;
        alpha_lo = alpha_prev;
        alpha_j = (alpha_lo + alpha_hi)/2;
        while abs(alpha_hi - alpha_lo) > 1e-6
            if feval(p,x+alpha_j*d,0) > feval(p,x,0)+i.c1ls*alpha_j*feval(p,x,1)'*d || (feval(p,x+alpha_j*d,0) >= feval(p,x+alpha_lo*d,0))
                alpha_hi = alpha_j;
            else
                phi = feval(p,x+alpha_j*d,1)'*d;
                if abs(phi) <= -i.c2ls*(feval(p,x,1)'*d)
                    alpha = alpha_j;
                    return;
                end
                if phi*(alpha_hi-alpha_lo) >= 0
                    alpha_hi = alpha_lo;
                end
                alpha_lo = alpha_j;
            end
            alpha_j = (alpha_lo+alpha_hi)/2;
        end
end
    