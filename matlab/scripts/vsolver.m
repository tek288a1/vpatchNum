function [a_term, rho_term, beta_term] = ...        
  vsolver(rho_init, rho_term, nr_step, which_beta, a_init, nr_term)
%VSOLVER [rho_term, beta_term, a_term] = 
%      vsolver(rho_init, rho_term, nr_step, which_beta, a_init, nr_term)
%
%  Newton's iteration method (continuation method) to solve for the
%  translational velocity  u  and the mapping coefficients of 2-D
%  vortex patch problem. rho is the continuation parameter.
%
%  INPUTS: 
%    rho_init = initial value of rho
%    rho_term = terminal value of rho; it can be smaller than rho_init.
%    nr_step = number of steps from rho_init to rho_term; nr_step + 1 
%              equispaced rho values including the initial and terminal 
%              ones are examined
%    which_beta = 0 : beta = 0
%               = 1 : beta = (1 - sqrt(1-rho^2))/rho
%               = 2 : beta = (1 - 2*sqrt(1-rho^2))/rho
%    a_init = initial iterate
%    nr_term = total number of unknowns (1 velocity + n coefficients); it
%              may be larger than the length of a_init.
%
%  NOTES: including logging using diary function
%
% Written by Tae Eun Kim, 12/11/2017

% stopping criteria
    tol_res = 1e-10;
    tol_rcnd = 1e-16;
    tol_dsol = 1e-10;
    max_iter = 10;
    rho_vec = linspace(rho_init, rho_term, nr_step+1);
    I = eye(nr_term);
    h = 1e-8;                           % Jacobian stepsize
    a0 = zeros(nr_term, 1);
    a0(1:length(a_init)) = a_init;
    for rho = rho_vec
        if which_beta == 0
            beta = 0;
        elseif which_beta == 1
            beta = (1-sqrt(1-rho^2))/rho;
        elseif which_beta == 2
            beta = (1-2*sqrt(1-rho^2))/rho;
        end    
        nr_iter = 0;
        a = a0;
        res = resval(rho, beta, a);
        fprintf('\n%7s %6.4f,', 'rho =', rho);
        fprintf('%9s %6.4f,', 'beta =', beta); 
        fprintf('%6s %5d\n', 'n =', nr_term);
        fprintf('%45s \n', repmat('-', 1, 43));
        fprintf('%6s %12s %12s %12s \n',...
                'iter', '|res|', '|dsol|', 'rcond');
        while true
            J = zeros(nr_term);         % Jacobian calculation
            for j = 1:nr_term
                res_pert = resval(rho, beta, a+h*I(:,j));
                J(:,j) = (res_pert - res)/h;
            end
            anew = a - J\res;           % iterative formula
            resnew = resval(rho, beta, anew);
            norm_res = norm(resnew, Inf);
            norm_dsol = norm(a - anew, Inf);
            rcnd = rcond(J);
            nr_iter = nr_iter + 1;
            fprintf('%6d %12.4e %12.4e %12.4e \n',...
                    nr_iter, norm_res, norm_dsol, rcnd);
            % stopping block
            if (norm_res <= tol_res) && (norm_dsol <= tol_dsol)
                a0 = anew;
                fprintf('%45s \n\n',...
                 '  ******** Solution has been found. *********');
                break;
            elseif nr_iter == max_iter
                fprintf('%45s \n\n',...
                 '  ********** No solution is found. **********');
                return;
            end
            a = anew;                   % updating iterates
            res = resnew;
        end
    end
    beta_term = beta;                   
    a_term = a0;
end
