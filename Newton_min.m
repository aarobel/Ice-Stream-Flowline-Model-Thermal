function [xout,error_flag,iter] = Newton_min(f,df,d2f,x,parameters,srchparams)
%Newton_min: Multidimensional Newton's method for function minimization,
%using line search step length determination
%
%   Input arguments:
%   f:      function handle of function to be minimized
%   df:     function handle of gradient of f
%   d2f:    function handle of Hessian of f
%   x:      initial guess; must have correct dimensions for an input
%           argument of f, df, d2f
%   parameters: parameter structure to be passed to f, df, d2f
%   srchparams: optional structure containing search parameters for Newton
%               iteration. This can contain the fields:
%               itmax:  maximum number of iterations
%               toldelta:   error bound on x; last Newton step must be <=
%                           0.5*toldelta
%               tolgrad:    error bound on df = 0; 0.5*norm(df,2) <=
%                           tolgrad to terminate iteration
%               tolsrch: if gradient is reduced by less than
%                           (1-tollinsrch), determine step length by line search
%                        0 < tolsrch < 1
%               verbose: if 0, runs quietly, if 1, outputs some
%                           information, if 2, lots of information
%
%   The algorithm used is Newton's method as a default, switching to a line
%   search using Brent's method if the gradient is not reduced
%   sufficiently (by at least a factor (1-tolsrch)) by a full Newton step

%NOTE: Check against Numerical Recipes for more efficient version; see
%below for comments regarding line search.

    %set up Newton iteration parameters
    if nargin == 6
        if isfield(srchparams,'itmax')
            itmax = srchparams.itmax;
        else
            itmax = 25;
        end
        if isfield(srchparams,'toldelta')
            toldelta = srchparams.toldelta;
        else
            toldelta = sqrt(eps);
        end
        if isfield(srchparams,'tolgrad')
            tolgrad = srchparams.tolgrad;
        else
            tolgrad = 100*eps;
        end
        if isfield(srchparams,'tolsrch')
            tolsrch = srchparams.tolsrch;
        else
            tolsrch = 5e-2;
        end
        if isfield(srchparams,'verbose')
            verbose = srchparams.verbose;
        else
            verbose = 0;
        end
    else
        itmax = 25;
        toldelta = sqrt(eps);
        tolgrad = 100*eps;
        tolsrch = 5e-2;
        verbose = 0;
    end

    %set parameters for line search
    paramsBrent.itmax = 100;
    paramsBrent.toldelta = toldelta;
    paramsBrent.bracketsize = 0.5;
    paramsBrent.dilation = 2;
    paramsBrent.verbose = 0;
    
    %set error flag to false if necessary
    if nargout == 3, error_flag = false; end

    %set up first Newton iteration step
    DF = df(x,parameters);
    DF_abs = tolgrad + eps;
    DF_abs_old = 0.5*norm(DF,2);
    if verbose, disp(strcat('Initial gradient =',num2str(DF_abs_old))), end
    Dx = toldelta + eps;
    iter = 1;
    while iter < itmax && (DF_abs >= tolgrad || norm(Dx) >= toldelta)
        D2F = d2f(x,parameters);
        Dx = -D2F\DF;
        x_temp = x + Dx;
        DF = df(x_temp,parameters);
        DF_abs = 0.5*norm(DF,2);
        if DF_abs > (1-tolsrch)*DF_abs_old
            %switch to line search
            if verbose, disp('Switching to line search'), end
            lambda = Brent_linsrch(f,df,x,Dx,parameters,paramsBrent);
                %NOTE: It may be overkill to find lambda to high accuracy
                %here; check numerical recipes and Dennis and Schnabel for 
                %better alternatives (Dennis and Schnabel have method
                %similar to Brent's but based purely on changes in gradient
                %rather than a tolerance on dlambda; should be possible to
                %amend Brent to that effect. Advantage of Brent of Dennis
                %and Schnabel is that Brent does not always use parabolic /
                %cubic interpolation, want to maintain that advantage.
                %Check numerical recipes for alternative.
                %NOTE: May also be possible to get Brent_linsrch to output
                %DF_abs. Consider re-coding.
            if verbose, disp(strcat('Step size =',num2str(norm(lambda*Dx,2)))), end
            x = x + lambda*Dx;
            DF = df(x,parameters);
            DF_abs = 0.5*norm(DF,2);
            if verbose, disp(strcat('Updated gradient =',num2str(DF_abs))), end
        else
            x = x_temp;
            if verbose, disp(strcat('Updated gradient =',num2str(DF_abs))), disp(strcat('Step size =',num2str(norm(Dx,2)))),  end
        end
%         if verbose, disp(strcat('Mean of solution =',num2str(mean(x)))), end
        if verbose, disp(strcat('Iteration Number =',num2str(iter))), end
        DF_abs_old = DF_abs;
        iter = iter + 1;
        if verbose == 2, disp('Current iterate ='), disp(x), disp('Current value of objective function:'), f(x,parameters), end
    end
    if iter >= itmax, warning('Newton minimization did not converge'), if nargout == 2, error_flag = true; end, end
    xout = x;
end

