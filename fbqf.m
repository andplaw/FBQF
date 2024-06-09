function [t,y] = fbqf(al,rhs,rhsp,t0,tf,y0,ht,b,k,stcl_max,tol,iterMax)

    % It's fine! Trust me! ------------------------------------------------
    id = 'MATLAB:nearlySingularMatrix';
    warning('off',id)
    
    id = 'MATLAB:SingularMatrix';
    warning('off',id)
    %----------------------------------------------------------------------

    % Check inputs --------------------------------------------------------
    if nargin <12
        iterMax = 5e3;
        if nargin < 11
            tol = 1e-14;
            if nargin < 10
                stcl_max = 7;
                if nargin < 9
                    k = 15;
                    if nargin < 8
                        b = 0.2;
                    end
                end
            end
        end
    end
    %----------------------------------------------------------------------
    %|
    %|
    % Setup ---------------------------------------------------------------

    %%% Exponential Change of Variables
    e2t = @(x) x.*exp(-b./x);
    de2t = @(x) exp(-b./x).*(b./(x+(x==0))+1);
    t2e = @(t) b./lambertw(b./t);

    %%% Monomial Change of Variables
    m2t = @(x) x.^k; 
    dm2t = @(x) k*x.^(k-1);
    t2m = @(t) t.^(1/k);

    %%% Miscellaneous
    N = round(tf/ht);
    y = zeros(1,N);
    
    %%% Initialize exponential node set
    e0 = t2e(t0);
    ef = t2e(tf);
    e_nodes = linspace(e0,ef,N);
    t = e2t(e_nodes);
    he = e_nodes(2) - e_nodes(1);

    %%% Identify transition point
    Te = e_nodes(stcl_max+2);   %transition point in exponential-time
    Tt = e2t(Te);               %transition point in untransformed-time
    Tm = t2m(Tt);               %transition point in monomial-time
    
    %%% Initialize monomial node set
    d = ceil(Tm/(N*(Tm - t2m(e2t(Te-he)))));    %Scaling factor
    Nm = d*N;                                   %number of monomial nodes
    m0 = t2m(t0);                               %starting point in mono-t
    m_nodes = linspace(m0,Tm,Nm);               %mono-t nodes up to transition
    hm = m_nodes(2) - m_nodes(1);               %mono-t node spacing
    ym = zeros(1,double(Nm));

    %%% Unit Tests
    unitTest_stepsize(m2t,m_nodes);
    unitTest_transitionPt(t2e,m2t,m_nodes,he)
    
    %----------------------------------------------------------------------
    %|
    %|
    % Solver --------------------------------------------------------------

    %================== Nonlinear Optimization Routine ===================%

    %%% Setup for nonlinear optimization
    polyO = 1;                              % Polynomial order to fit
    kickstarter_r_ep = 7;                   % Number of points to kickstart
    unitTest_underflow(kickstarter_r_ep,m_nodes,m2t,al);

    %%% Initialize domain w/ nonlinear optimization routine
    y_ks = kickstarter(m_nodes(1:kickstarter_r_ep),y0,polyO,k,al,rhs);
    ym(1:kickstarter_r_ep) = y_ks;


    %============================ FBQF (mono) ============================%

    %%% Setup base & eval point stencils and weights
    stcl_bp = 1:stcl_max;                   % Base point stencil
    m_bp = m_nodes(stcl_bp);                % Base point stencil nodes
    cwbp = EC_weights(-1,hm,m_bp);          % Base point stencil weights
    cwep = fliplr(EC_weights(al,hm,m_bp));  % Eval point stencil weights

    %%% Step FBQF forward in (monomial) time
    for i = (kickstarter_r_ep):Nm-1

        stcl_ep = i-stcl_max+2:i+1;         % Eval point stencil nodes

        [S,wep0] = fbqfSetup(m_nodes(1:i+1),ym(1:i),...
                             cwbp,cwep,...
                             stcl_bp,stcl_ep,...
                             al,hm,...
                             @(x)m2t(x),@(x)dm2t(x));

        ym(i+1) = fbqfNewtons(S,wep0,...
                            m2t(m_nodes(i+1)),ym(i),...
                            rhs,rhsp,...
                            al,...
                            tol,iterMax);

        unitTest_NaN(ym(i+1),"ym(i+1)")
        unitTest_Imag(ym(i+1),"ym(i+1)")
        
    end
    

    %=========================== Interpolation ===========================%
    
    %%% Interpolate
    p = griddedInterpolant(m2t(m_nodes),ym(1:Nm));
    y_interp = p(e2t(e_nodes(1:(stcl_max))));

    %%% A bit of care is necessary when negative y values upset the rhs of
    %%% the FODE. The interpolation may result in some function values 
    %%% becoming negative for very small step sizes.
    if all(y_ks>=0)&&any(y_interp<0)
        y_interp = abs(y_interp);
    end

    %%% Commit interpolant points
    y(1:stcl_max) = y_interp; 
    

    
    %============================= FBQF (exp) ============================%

    %%% Setup base & eval point stencils and weights
    % Remains the same                      % Base point stencil
    % Remains the same                      % Base point stencil nodes
    % Remains the same                      % Base point stencil weights
    % Remains the same                      % Eval point stencil weights
    
    %%% Step FBQF forward in (exponential) time
    for i = (stcl_max):N-1

        stcl_ep = i-stcl_max+2:i+1;         % Eval point stencil nodes

        [S,wep0] = fbqfSetup(e_nodes(1:i+1),y(1:i),...
                            cwbp,cwep,...
                            stcl_bp,stcl_ep,...
                            al,he,...
                            @(x)e2t(x),@(x)de2t(x));

        y(i+1) = fbqfNewtons(S,wep0,...
                            e2t(e_nodes(i+1)),y(i),...
                            rhs,rhsp,...
                            al,...
                            tol,iterMax);
        
    end

end


%=========================================================================
%
% HELPER FUNCTIONS
%
%=========================================================================

function yv = kickstarter(xv,y0,p,m,al,rhs)

    arguments
        xv      % Points in transformed domain
        y0      % Initial condition
        p       % Order of polynomial approximation
        m       % Order of monomial change of variables
        al      % Order of fractional derivative
        rhs     % Right hand side of the FODE
    end
    
    % This function builds the messy entries of the A matrix
    % The solution to the exact Caputo derivative of a monomial
    helper = @(k,r) gammaDiv((1+k)/m,1-al+(1+k)/m)./m.*r.^(1+k-m*al);

    % Function for transformation of indep var
    x2t = @(x) x.^m; 

    A = zeros([numel(xv)-1,p]);
    for k = 0:p
        for ix = 2:numel(xv)
            A(ix-1,k+1) = helper(k,xv(ix));
        end
    end
    
    unitTest_overflow(A)

    yVand = zeros(length(xv)-1,double(p+1));
    for k = 0:p
        yVand(:,k+1) = xv(2:end)'.^(k+1)./(k+1);
    end

    diffun = @(a) A*a - rhs(x2t(xv(2:end))',yVand*a+y0.*ones(length(xv)-1,1));
    if log10(max(max(A)))>150 %overflow will be reached
        scale = 2^((log(min(min(A)))+log(max(max(A))))/(log(2)*2));
        diffexpr = @(a) diffun(a)/scale;
    else 
        diffexpr = @(a) diffun(a);
    end

    xInit.a = zeros(1,double(p+1));

    options = optimoptions("lsqnonlin");
    options = optimoptions(options,"FunctionTolerance",1e-14,...
                                       "ConstraintTolerance",1e-14,...
                                       "StepTolerance",1e-14,...
                                       "OptimalityTolerance",1e-14,...
                                       "MaxIterations",2*options.MaxIterations, ...
                                       "Display","off"); %iter-detailed
    [aCoeffs,~,~,~,~] = lsqnonlin(diffexpr ,...
                              xInit.a',...
                              [],...
                              [],...
                              options);
    
    yCoeffs = polyint(flipud(aCoeffs)',y0);
    yv = polyval(yCoeffs,xv(1:end));

    function prod = gammaDiv(x1,x2)
    
        prod = 1;
        for i = 1:min(floor(x1)-1*(x1==floor(x1)),floor(x2)-1*(x2==floor(x2)))
            
            factor = (x1-i)/(x2-i);
            prod = factor*prod;
        end
    
        if isempty(i)
            i=0;
        end
    
        prod = gamma(x1-i)/gamma(x2-i)*prod;
    end

end

function w = EC_weights(al,h,z) 
    %   Calculate weights in a stencil matching the Taylor expansion 
    %   ζ(α+1)/0!+ζ(α)(hξ)/1!+ζ(α-1)(hξ)²/2!+ζ(α-2)(hξ)³/3!+ζ(α-3)(hξ)⁴/4!+...    
    % Input parameters 
    %  al   The parameter denoted α above 
    %   h   Step size in TR-type approximation (with first term omitted) 
    %   z   Stencil node locations (array or vector) 
    % Output parameter 
    %   w   Weights at locations given by z  
    
    dm = size(z);                    % Dimensions of z; to use in reshape for w  
    z = z(:);  sq = (0:length(z)-1); % Make z a column vector 
    rhs = zeta(al+1-sq).*h.^sq;      % Get RHS vector for weight calculation 
    A = rot90(vander(z));            % vandermonde coefficient matrix 
    w = A\rhs';                      % Solve linear system 
    w = reshape(w,dm);               % Arrange weights in same way used for z 
    w = double(w);
end 

function [S,wep0] = fbqfSetup(x,y,cwbp,cwep,stcl_bp,stcl_ep,al,h,eta,eta_p)

    arguments
        x               %indep function values - with transformation
        y               %dependent function values
        cwbp            %coeff of g_ns from base correction stencil
        cwep            %coeff of g_ns from end  correction stencil
        stcl_bp         %base correction stencil
        stcl_ep         %end  correction stencil
        al              %order of fractional derivative
        h               %step size ante-transition
        eta             %variable change
        eta_p           %derivative of variable change
    end
    
    N = length(x);
    ep = x(N); %end point
    
    %%% transformed variable after CoV and IBP
    g = @(s) eta_p(s)./(eta(ep)-eta(s)).^(al+1); %factor from integration by parts
    gv = g(x(1:end-1)); 
    gv = [gv,eta_p(x(end))^(-al)];    % lim_{s->lambda} g1(s) = 1/eta_p(lambda)

    
    %%% Left End Point - bp
    ywbp = y(stcl_bp);
    gvbp = gv(stcl_bp);
    
    %%% TR
    TR1 = h*sum(gv(2:end-1).*y(2:end));    % Evaluate the TR-type sum 
    
    %%% Right End Point - ep
    % The final weight is multiplied by the transformation terms so that 
    % the calculations don't need to be redone outside of this function
    wep0 = h^(-al)*cwep(end)*gv(end); 
    ywep = y(stcl_ep(1:end-1)); %+2 = +1 for the unknown and +1 for counting
    gvep = gv(stcl_ep(1:end-1)).*(ep-x(stcl_ep(1:end-1))).^(al+1);
    
    
    %%% Output
    % full result without y_{n+1} term in end point stencil (since it's unknown)
    
    S = -y(1)/(eta(ep)^al)...
        - al*(TR1...
              - h*(cwbp.*ywbp)*gvbp'...
              - h^(-al)*(cwep(1:end-1).*ywep)*gvep'); 

    unitTest_Imag(S,'S')
    unitTest_NaN(S,'S')
end


function [ynp1] = fbqfNewtons(S,wep0,ep,yold,rhs,rhsp,al,tol,iterMax)

    arguments
        S               %previous function values in BDF method
        wep0            %coefficient of unknown y value
        ep              %value of end point in rhs input variable
        yold            %dependent function value y(n)
        rhs             %right hand side of FDE
        rhsp            %right hand side of FDE, differentiated
        al              %order of fractional derivative
        tol = 1e-16     %tolerance dictating when iterations end
        iterMax = 5e3   %max iterations if tolerance is not reached
    end

    %%% Newton's Method
    u = @(g) (S + al*wep0*g)/gamma(1-al) - rhs(ep,g);
    up = @(g) al*wep0/gamma(1-al) - rhsp(ep,g);
    
    ynp1 = yold - u(yold)/up(yold);
    
    changed = true;
    j = 0;

    while abs(ynp1-yold)>tol && changed
        yoldold = yold;
        yold = ynp1;
        ynp1 = yold - u(yold)/up(yold);
        j = j+1;
        changed = logical(abs(ynp1-yold)<abs(yold-yoldold));

        if j==iterMax
            msgStr = "Newtons method reached the max iterations before ";
            msgStr = msgStr + "reaching the desired tolerance.\n";
            warning(msgStr)
            break
        end
    end

    unitTest_Imag(ynp1,"y(i+1)")
end

%=========================================================================
%
% UNIT TESTS
%
%=========================================================================

function unitTest_underflow(kickstarter_r_ep,m_nodes,m2t,al)

    % Integration by parts term becomes infinite
    if (m2t(m_nodes(kickstarter_r_ep+1))-m2t(m_nodes(kickstarter_r_ep)))^(1+al)==0
        msgStr = "The transformation of the monomial nodes into t nodes ";
        msgStr = msgStr + "then taken to the (al+1) power is too ";
        msgStr = msgStr + "small. Decrease monomial power, alpha or increase";
        msgStr = msgStr + "step size to resolve this issue.";
        error(msgStr)
    end

    
end

function unitTest_stepsize(m2t,m_nodes)

    if find(m2t(m_nodes)>0,1) > 5
        
        msgStr = "The step size is too small in the monomial domain.\n";
        msgStr = msgStr + "Please increase the step size.\n";
        %{
            This error occurs when the conversion from the m nodes into the
            t variable produces too many t values below underflow. This is
            the result of the step size being too small.
        %}
        % msgStr = msgStr + "You can also try decreaseing b in the ";
        % msgStr = msgStr + "exponential change of variables.";
        error(sprintf(msgStr))
    end
end

function unitTest_overflow(A)

    % optimization terms in jacobian multiplication become infinite
    if log10(max(max(abs(A))))>150
        scale = 2^((log(min(min(A)))+log(max(max(A))))/(log(2)*2));
        A = A./scale;
        if isinf(max(max(abs(A)))^2)
            msgStr = "The transformation of the monomial nodes into t nodes ";
            msgStr = msgStr + "makes them too ";
            msgStr = msgStr + "small. Decrease monomial power, alpha or increase";
            msgStr = msgStr + "step size to resolve this issue.";
            error(msgStr)
        end
    end
end

function unitTest_NaN(query,string)

    if isnan(query)
        msgStr = string + " should not be NaN";
        error(msgStr)
    end
end

function unitTest_Imag(query,string)

    if imag(query)~=0
        msgStr = string + " should not be imaginary";
        error(msgStr)
    end
end

function unitTest_transitionPt(t2e,m2t,m_nodes,he)

    if (diff(t2e(m2t(m_nodes(end-1:end))))/he-1)>1e-1
        msgStr = "The step size at the transition point isn't ";
        msgStr = msgStr + "the same in both the monomial CoV ";
        msgStr = msgStr + "and the exp CoV.";
        error(msgStr)
    end
end