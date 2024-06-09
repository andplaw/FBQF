%% Preamble

close all
clear all

set(groot,'defaultLineLineWidth',1.5)
set(groot,'defaultAxesFontSize',12)
set(groot,'DefaultAxesFontWeight','bold');
set(groot,'defaultAxesTitleFontSizeMultiplier',1.33)
set(groot,'defaultAxesTitleFontWeight','bold')



bool_save = false; 

%% Parameters

%%% Order of fractional derivative
% al = 0.5; %seems to converge more poorly for low alpha
als = 0.5; % [0.01, 0.1,0.25,0.5,0.75,0.9,0.99]; %+(pi-3)/100; %0.1+(pi-3)/100; %pi-3+0.4; %0.5; %[0.01, 0.1,0.25,0.5,0.75,0.9,0.99]; %

%%% time interval bounds
t0 = 0;
tf = 1;

%%% Number of steps across interval
Ns = round(logspace(log10(2^8),log10(2^12),15))+1; %2.^(6:12)+1; %2.^(4:12)+1; %[64 128 256 512 1024 1536 2048 3072 4096 6144 8192] + 1; %1+round(logspace(1.3,3,10)); %
hs = tf./(Ns);

%%% convergence tolerance for Newton's method
tol = 10^(-16);

%%% reference line for convergence plots
h_ref = 6; %examples{ie,4};


%% Examples
% some examples are true for any alpha, 0<al<1. Those have some explicit
% dependence on alpha, accordingly
% some examples are true only for alpha, al==1/2. These do not have
% explicit alpha dependence

y0 = 1;
y02 = 0;
lambda = 1; 
f = @(t,al) (gamma(6)/gamma(6-al))*t.^(5-al)-3*(gamma(5)/gamma(5-al))*t.^(4-al)+...
         2*(gamma(4)/gamma(4-al))*t.^(3-al)+(t.^5-3*t.^4+2*t.^3).^2;
f2 = @(t,al) (gamma(2*al+1)/gamma(al+1)).*t.^al - 2/gamma(3-al).*t.^(2-al) + (t.^(2*al)-t.^2).^4;
f3 = @(t,al) 40320/gamma(9-al).*t.^(8-al) - 3*gamma(5+al/2)/gamma(5-al/2).*t.^(4-al/2)+...
          9/4*gamma(al+1) + (3/2*t.^(al/2)-t.^4).^3;
f4 = @(t,al) t.^4 - 1/2.*t.^3 - 3/gamma(4-al).*t.^(3-al) + 24/gamma(5-al).*t.^(4-al);
f5 = @(t,al) 40320/gamma(9-al).*t.^(8-al) - 3*gamma(5+al/2)/gamma(5-al/2).*t.^(4-al/2)+...
          9/4*gamma(al+1) ;

nu1 = sqrt(2)*10^(-2); %sqrt(2)/140
nu2 = 1-nu1; %sqrt(2)/2+0.292;
nu3 = 0.5; %1*sqrt(2)/6*10^(-1);
v1 = 0;
v2 = sqrt(2)/2+0.29;
     
examples_al = @(al) {...
         @(t) t.^5 - 3*t.^4 + 2*t.^3, @(t,y) f(t,al) - y.^2, @(t,y) -2.*y, 7, "$y(t) = t^5 - 3t^4 + 2t^3$"; %1
         @(t) y0.*mlf(al,1,-lambda.*t.^al,mlf_p)-1, @(t,y) -lambda.*(y+1), @(t,y) -lambda, 6, sprintf("$y(t) = E_{\\alpha}(%g t^{\\alpha})$",-lambda); %2
         @(t) t.^(2*al)-t.^2, @(t,y) f2(t,al) - y.^4, @(t,y) -4.*y.^3, 9, "$y(t) = t^{2\alpha} - t^2$"; %3
         @(t) t.^2, @(t,y) t.^2 - y + 2/gamma(3-al).*t.^(2-al), @(t,y) -1, 9, "$y(t) = t^2$"; %4
         @(t) t.^8 - 3*t.^(4+al/2)+9/4*t.^al, @(t,y) f3(t,al)-y.^(3/2), @(t,y) -3/2*y.^(1/2), 7,"$y(t) = t^8 - 3t^{4+\alpha/2}+\frac{9}{4}t^{\alpha}$"; %5
         @(t) 1+(y02-1).*mlf(al,1,-lambda.*t.^al,mlf_p)-y02, @(t,y) lambda.*(1-(y-y02)), @(t,x) -lambda, 7, sprintf("$y(t) = 1+ %g E_{\\alpha}(%g t^{\\alpha})$",y02-1,-lambda); %6
         @(t) (t.^(al+1) - t.^al)./gamma(al+1), @(t,y) (1+al).*t-1, @(t,y) 0, 6, "$y(t) = \frac{t^{\alpha+1}-t^\alpha}{\Gamma(\alpha+1)}$"; %7
         @(t) 1-mlf(al,1,t.^al,mlf_p)+t.^(1+al).*mlf(al,2+al,t.^al,mlf_p), @(t,y) t+(y-1), @(t,y) 1, 8, "$y(t) = t^{\alpha+1} E_{\alpha,\alpha+2}\left(t^\alpha\right)- E_{\alpha}\left(t^\alpha\right)+1$"; %8
         @(t) t.^4 - 1/2.*t.^3, @(t,y) f4(t,al) - y, @(t,x) -1, 6, "$y(t) = t^4 - \frac{1}{2} t^3$"; %9 [Al-Rabatah 2012]
         @(t) 1 - mlf(al,1,-t.^al,mlf_p) + t.^(1+al).*mlf(al,2+al,-t.^al,mlf_p), @(t,y) t-(y-1), @(t,y) -1, 6, "$y(t) = t^{\alpha+1} E_{\alpha,\alpha+2}\left(-t^\alpha\right)- E_{\alpha}\left(-t^\alpha\right)+1$"; %10
         @(t) t.^nu1, @(t,y) (t.^(nu1-al).*gamma(1+nu1))./gamma(nu1-al+1),@(t,y) 0,6, sprintf("$y(t) = t^\\nu, \\nu=%.4f$",nu1); %11
         @(t) t.^v2 + t.^(v1*al), @(t,y) t.^(v2-al)*gamma(1+v2)/gamma(1+v2-al)+t.^((v1-1)*al)*gamma(1+v1*al)/gamma(1+(v1-1)*al),@(t,y) 0,6, sprintf("$y(t) = t^{%g}+t^{v_1 \\alpha}, v_1=%g$",v2,v1) %12
         @(t) t.^8 - 3*t.^(4+al/2)+9/4*t.^al, @(t,y) f5(t,al), @(t,y) 0, 7,"$y(t) = t^8 - 3t^{4+\alpha/2}+\frac{9}{4}t^{\alpha}$"; %13
         @(t) t.^nu2, @(t,y) (t.^(nu2-al).*gamma(1+nu2))./gamma(nu2-al+1),@(t,y) 0,6, sprintf("$y(t) = t^\\nu, \\nu=%.4f$",nu2); %14
         @(t) t.^al, @(t,y) (gamma(1+al))./gamma(1),@(t,y) 0,6, "$y(t) = t^\alpha$"; %15
         @(t) t.^nu3, @(t,y) (t.^(nu3-al).*gamma(1+nu3))./gamma(nu3-al+1),@(t,y) 0,6, sprintf("$y(t) = t^\\nu, \\nu=%.4f$",nu3)}; %16


         % @(t) -4.*t.^(3/2)./(3*sqrt(pi)), @(t,y) -t, @(t,y) 0, 7, "$y(t) = \frac{4t^{3/2}}{3\sqrt{\pi}}$"; %2
         % @(t) t.^(1/2), @(t,y) gamma(3/2), @(t,y) 0, 6, "$y(t) = t^{1/2}$"; %3
 
%% FDE Solver - Backward Differentiation


%%% Initialization
errs = zeros(size(examples_al,1),size(als,2),size(Ns,2));
errs_end = zeros(size(examples_al,1),size(als,2),size(Ns,2));
tocs = zeros(size(examples_al,1),size(als,2),size(Ns,2));

for ie = 1 %[1,2,11,14,16] %1:size(examples_al(NaN),1)

    for ial = 1:size(als,2)
        
        al = als(ial);
        examples = examples_al(al);
        y_exact = examples{ie,1};
        rhs = examples{ie,2};
        rhsp = examples{ie,3};

        for iN = 1:size(Ns,2)

            %initialize steps/step-size
            N = Ns(iN);
            ht = hs(iN);

            tic

            [t,y] = fbqf(al,rhs,rhsp,t0,tf,y0,ht);
            
            %%% Record elapsed time
            tocs(ie,ial,iN) = toc;

            %%% Quantify the error
            errs(ie,ial,iN) = max(abs(y_exact(t(1:end))-y(1:end)));
            errs_end(ie,ial,iN) = abs(y_exact(t(end))-y(end));
        end
        
        
    end

    %%% Legend entries
    if iscell(als)
        for iik = 1:size(als,1)
            ks_cell{iik} = als(iik,:);
        end
    else
        ks_cell = num2cell(als);
    end
    legend_entries = cellfun(@(x) "al="+num2str(x),ks_cell,UniformOutput=false);
    legend_entries{end+1} = sprintf('h^{%.2f}',h_ref);


    %%% Plot max abs error and absolute error at t=tf
    errs_xs = {errs,errs_end};
    errs_xs_str = {"Maximum Absolute Error",sprintf("Absolute Error at t=%g",tf)};
    errs_xs_abbv_str = {"max_err","end_err"};

    for ierrs_x = 1:numel(errs_xs)

        errs_x = errs_xs{ierrs_x};
        ev_N = reshape(errs_x(ie,ial,:),1,size(Ns,2));
        errs_xs_str{ierrs_x}
        
        ev = reshape(errs_x(ie,:,:),size(als,2),size(Ns,2));
        EOC = log(ev(:,2:end)./ev(:,1:end-1))...
              ./log(hs(2:end)./hs(1:end-1))
        
        figure
        loglog(hs,ev')
        hold on
        loglog(hs,hs.^(h_ref)*10^(round(log(ev(1,1))/log(10)-log(hs(1).^(h_ref))/log(10))),'--k')
        
        ylim([1e-16 1e0])
        xlabel('Step Size ($h_\theta$)',"interpreter","latex")
        ylabel(errs_xs_str{ierrs_x},"interpreter","latex")
        set(gca,'xdir','reverse')%,'xscale','log','zscale','log')
        % axis tight
        axis square
        title(examples{ie,5},"interpreter","latex") %+" $\alpha ="+num2str(al)+"$ $\beta = "+num2str(b)+"$ $k = "+num2str(k)+"$"
        legend(legend_entries,'Location','southwest')
        drawnow
        
        if bool_save
            ax = gca();
            filename = "Figures/FDENLS_Example_"+num2str(ie)+"_b_"+num2str(b)+"_k_"+num2str(k)+"_"+errs_xs_abbv_str{ierrs_x}; %+"_al_"+num2str(al)
            exportgraphics(ax,filename+".eps",'BackgroundColor','none')
            exportgraphics(ax,filename+".pdf",'BackgroundColor','none')
            close
        end
    end

end


% [1] E.A.Rakhmanov,Bounds for polynomials with a unit discrete norm, Ann. Math. 165 (2007), 55–88.
%       on convergence/stability tradeoff of least-squares polynomial
%       fitting
% [2] J.P.Boyd and F. Xu, Divergence (Runge phenomenon) for least-squares polynomial approximation on an equispaced grid and Mock-Chebyshev subset interpolation, Appl. Math. Comp. 210 (2007), 158–168.
%       numerical experiments on convergence/stability tradeoff  of
%       least-squares polynomial fitting

