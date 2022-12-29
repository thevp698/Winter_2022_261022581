function x= optsolver(p,x,algorithm,i)
i.maxiter = 1e+3;   %number of iterations
i.opttol = 1e-6;    %tolerance of the gradient 
i.c1ls = 0.30;      %C1 for wolfe condition
i.c2ls = 0.70;      %C2 for wolfe condition
i.c1tr = 0.10;      %C1 in the trust region
i.c2tr = 0.90;      %C2 in the trust region
i.cgopttol = 1e-6;  %used while solving the trust region subproblem AX = B
i.cgmaxiter = 1e+3; %USed while solving the trust region subproblem AX = B
i.sr1updatetol = 1e-6;%SR1 update tolerance
i.bfgsdamptol = 1e-6; %BFGS update tolerance
f = feval(p,x,0);   %function value at the intial value of x
g = feval(p,x,1);   % gradient value at the intial value of x
H = feval(p,x,2);   % hessian value at the intial value of x
tr = 2;             % intializing the radius of the trust region
% I tried create simple program to find the optimum of the function by
% calculatin the direction and the stepsize and calculate the next iterate
% by x = x + ad. Our function reaches to optimum point when the gradient of
% the function is very small or it reaches the maximum iterations.
%I've made four functions first one is 'optsolver' in
%which I call three other functions 'search direction',
%'stepacceptance','hessianupdate'
%in the optsolver I print norms, direction, number of iterations and
%stepsize.
% in the 'searchdirection' fucntion we find the direction of the function
% by differnt algoritms such as steepestdescent, Newton, CG and BFGS.
% In the stepsize function we take the appropriate stepsize as per the
% alorithms such as backtrack, wolfe, trusrtregion, SR1trustregion.
% In the hessian update function I update hessian as per the algoritms are being used, for BFGS and SR1 we update our hessian by iterative methods other than those two algoritms I calculated hessian by feval(p,x,2).  
% Store output strings
out_line = '==================================================================';
out_data = '  k        f          ||g||       ||d||       g^T*d        alpha';
out_null =                                '----------  -----------  ----------';

% Print output header
fprintf('%s\n%s\n%s\n',out_line,out_data,out_line);

% Initialize iteration counter
k = 0;

% Evaluate gradient norm
norms.g = norm(g);

% Store initial gradient norm
norms.g0 = norms.g;
tic
% Iteration loop
while 1
  
  % Print iterate information
  fprintf('%4d  %+.4e  %.4e  ',k,f,norms.g);
  
  % Check termination conditions
  if k > i.maxiter || norms.g <= i.opttol*max(norms.g0,1), break; end  
   
  % Evaluate search direction
  [d,D] = searchdirection(g,H,i,tr,algorithm); 
  
  % Evaluate norm of direction 
  norms.d = norm(d);
  
  fprintf('%.4e  ',norms.d);
  
  % Evaluate stpe acceptance
  [a,f,tr] = stepacceptance(p,x,f,d,D,i,H,tr,algorithm);
  
  % Evaluate directional derivative
  D = g'*d;
  
  % Print line search information
  fprintf('%+.4e  %.4e\n ',D,a);
  
  % Update iterate
   x= x+a*d;      
  
  %Hessian Update after the update of the iterate here I've defined one
  %function for the hessian where I update the hessian as per the algoritm
  %like if the algorithm is SR1 or BFGS then I'll use the iterative
  %otherwise I'll calculate hessian as H = feval(p,x,2)
  H = hessianupdate(H,p,x,a,d,g,algorithm,i);
  
  % Evalaute objective functiom
  f = feval(p,x,0);
  
  % Evaluate objective gradient
  g = feval(p,x,1);

  % Evaluate gradient norm
  norms.g = norm(g);
   
  % Increment iteration counter
  k = k + 1;

  %Store the value of k in an iteration to plot the graph
  iter = k; 

  %Store the function values in array
  fvalues(iter) = f;
  
  %Store the values of norms of gradient in an array
  normG(iter) = norms.g;

end
toc
% Print output footer
fprintf('%s\n%s\n',out_null,out_line);
% Convergence graph
% Creating a linear graph of norm of gradient VS number of iteration
figure
plot(normG)
title([algorithm ' '  'linear - linear axis plot'])
xlabel('iteration count')
ylabel('gradient norm')

%creating semilog graph
figure
semilogy(normG)
title([algorithm  ' ' 'linear-log axis plot'])
xlabel('iteration count')
ylabel('gradient norm')

%Creating a log-log graph0.3
figure
loglog(normG)
title([algorithm  '  ' 'log-log axis plot'])
xlabel('iteration count')
ylabel('gradient norm')

%function values VS number of iterations
figure
plot(fvalues)
title([algorithm])
xlabel('iteration count')
ylabel('function values')