function xout = Brent_linsrch(f,df,x0,v0,parameters,srchparams)
%function xout = Brent_linsrch(f,df,x0,v0,parameters,srchparams)
%Unconstrained pointwise minimization routine along a given search
%direction v0 from starting point x0 by Brent's
%method with gradient evaluation (see Numerical Recipes). The minimizer
%lies at x0 + xout*v0. A step length of one is assumed as an initial
%guess. Input format
%   f:  function handle for function to be minimized
%   df: function handle for gradient of function to be minimized;
%       the derivative with respect to step length x is then df*v0
%   x0: column vector giving starting point for line search
%   v0: search direction, same size vector as x0
%       f and df take arguments of size(x0)
%   parameters: parameter structure to be passed to f and df
%   srchparams: parameters for line search algorithm, optional
%   may contain fields
%       itmax: maximum number of iterations
%       toldelta: error bound on x
%       bracketsize: initial bracket size
%       dilation: dilation factor for brackets (> 1)
%       verbsoe: set to unity for verbose output

    %set up search parameters
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
        if isfield(srchparams,'bracketsize')
            bracketsize = srchparams.bracketsize;
        else
            bracketsize = .5;
        end
        if isfield(srchparams,'dilation')
            dilation = srchparams.dilation;
        else
            dilation = 2;
        end
        if isfield(srchparams,'verbose')
            verbose = srchparams.verbose;
        else verbose = 0;
        end
    else
        itmax = 25;
        toldelta = sqrt(eps);
        bracketsize = .5;
        dilation = 2;
        verbose = 0;
    end

    DF = df(x0+v0,parameters).'*v0;
    %line search assumes full step length as initial guess (a = 1); bracket
    %from there
    a = 1;
    while 1
        DF_old = DF;
        b = a-sign(DF)*bracketsize;
        DF = df(x0+b*v0,parameters).'*v0;
        if DF*DF_old > 0
            a = b; bracketsize = bracketsize*dilation;
        else
            if a > b, atemp = b; b = a; a = atemp; end
                if verbose, disp(strcat('Bracket using dilations at location ',num2str(ii),'a = ',num2str(a),'b = ',num2str(b))), end
                break
        end
    end

    %Now start Brent iteration
    iter = 1;
    v = (a+b)/2; w = v; x = v;          %during iteration, v, w = last two `best' guesses, x = current best guess
    Fx = f(x0+x*v0,parameters); Fv = Fx; Fw = Fx;
    DFx = df(x0+x*v0,parameters).'*v0; DFv = DFx; DFw = DFx;
    Dx_old = 0;
    while iter < itmax
        if abs(Dx_old) > toldelta      %step before last not too small, use parabolic interpolation with gradient
            if DFx ~= DFw, Dx1 = (w-x)*DFx/(DFx-DFw); else Dx1 = 2*(b-a); end
            if DFx ~= DFv, Dx2 = (v-x)*DFv/(DFx-DFv); else Dx2 = 2*(b-a); end
            x1 = x+Dx1; x2 = x+Dx2;     %check which updated estimates lie in bracket and which steps are in descent direction
            check1 = (DFx*Dx1 < 0 || (x1-a)*(b-x1) > 0);
            check2 = (DFx*Dx2 < 0 || (x2-a)*(b-x2) > 0);
            Dx_old2 = Dx_old; Dx_old = Dx;      %update previous step sizes
            if check1 || check2                 %pick smaller allowed step; if none are valid then do bisection instead
                if check1 && check2
                    Dx = -min(abs(Dx1),abs(Dx2))*sign(DFx); %NOTE: CHECK AGAINST NUMERICAL RECIPES
                elseif check1
                    Dx = Dx1;
                else
                    Dx = Dx2;
                end
                if abs(Dx) > 0.5*Dx_old2 ||   abs(x+Dx-a) < toldelta ||  abs(x+Dx-b) < toldelta     %only use parabolic interpolation if step size is less than half step size from two steps ago (why not last step?)
                    Dx = -toldelta*sign(DFx); Dx_old = Dx;     %do not allow parabolic interpolation result too close to bracket boundary, this does not update x and forces bisection on next step
                    if verbose == 2, disp('Brent parabolic interpolation unsuccessful, new iterate too close to bracket boundary'), end
                elseif verbose == 2, disp('Brent parabolic interpolation successful')
                end
            else
                Dx = -toldelta*sign(DFx); Dx_old = Dx;         %as above, bad parabolic interpolation forces bisection on next step
                                                                     %NOTE: CHECK AGAINST NUMERICAL RECIPES
                if verbose == 2, disp('Brent parabolic interpolation unsuccessful, new iterate too close to bracket boundary'), end
            end
        else
            if DFx > 0              %bisect interval in descent direction
                Dx_old = a-x;
            else
                Dx_old = b-x;
            end
                Dx = 0.5*Dx_old;
                if verbose == 2, disp('Bisect in descent direction'), end                    
        end
        if abs(Dx) > toldelta
            u=x+Dx;
            Fu = f(x0+u*v0,parameters);
            if verbose == 2, disp('Full step length taken'), end
        else
            u=x+toldelta*sign(Dx);     %enforce minimum step
            Fu = f(x0+u*v0,parameters);
            if verbose == 2, disp('Minimum step length taken'), end                   
            if Fu > Fx, 
                xout = x;
                if verbose == 2, disp('Minimum attained after minimum step is taken'), end
                break   %exit here if minimum step leads to increase in f 
            end
        end
        DFu = df(x0+u*v0,parameters).'*v0;
        if Fu <= Fx                     %Make u new best guess, update, v and w.
            if u > x                    %update bracketing interval
                a = x;
            else
                b = x;
            end
            v = w; Fv = Fw; DFv = DFw;
            w = x; Fw = Fx; DFw = DFx;
            x = u; Fx = Fu; DFx = DFu;
        else                            %Keep x as best guess, decrease bracketing interval
            if u < x                    %update bracketing interval
                a = u;
            else
                b = u;
            end
            if Fu <= Fw || w == x       %Now update one of v or w if this is warranted, i.e. if u is a better guess than v or w, or if at least one pair of v,w,x are equal
                v = w; Fv = Fw; DFv = DFw;
                w = u; Fw = Fu; DFw = DFu;
            elseif Fu <= Fv || v == x || v == w
                v = u; Fv = Fu; DFv = DFu;
            end
        end
        if abs(x-0.5*(a+b)) < 2*toldelta - 0.5*(b-a), 
            xout = x;
            if verbose == 2, disp('Minimum attained by shrinking bracket'), end
            break 
        end    %terminate if bracketing interval has become sufficiently small
        iter = iter + 1;
    end
    if iter >= itmax, warning(strcat('Too many iterations in Brent algorithm; a = ',num2str(a),' b = ',num2str(b),' x = ',num2str(x))), xout = x; end
end