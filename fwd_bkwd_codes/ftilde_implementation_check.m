% Total stress

% Create axes
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
%axes1.XScale='log';
%axes1.YScale='log';
title('$I_k$ for 3D multi-bead-spring-dashpot, $N=30$','Interpreter','latex','FontSize',20);
% xlabel('$t^*$','FontSize',30,'Interpreter','latex');
xlabel('Calculated from $f$','FontSize',30,'Interpreter','latex');
y=ylabel('Calculated from $\tilde{f}$','FontSize',36,'Interpreter','latex',...
    'Rotation',90);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
% xlim([0. 5.]);
% ylim([0. 1.]);
pbaspect([1. 1. 1.]);
format long;
% % grid on;


%number of springs
N=30;
%number of dimensions
ndim=3;
%Value of "k". The "M" value that is evaluated
k=8;
%value of "j". Connector vector with respect to which gradient is measured
j=7;

% rng(16032020);
% rng(400042020);
rng('shuffle');

%\omega=varphi/zeta
varphi=200;

p=(varphi/((2*varphi)+1))^2;

%creating initial configurations
Q=normrnd(0,1,[N,ndim]);
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  CALCULATING I_k using "f"            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_k_fromf = zeros(N,1);

for k=1:N
    I_k_fromf(k)=fwd_coeff_poly(k,L,p,N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  CALCULATING I_k using "ftilde"       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_k_fromftilde = ones(N,1);

for k=1:N
    for nu=1:(k-1)
        compare_array=gen_list(nu); %this involves repetitive calculation
                                    %room for speeding up
        minval=min(nu);
        I_k_fromftilde(k)=I_k_fromftilde(k)+(L(nu)*L(nu)*fwd_coeff_lsquared_nu_composite(nu,compare_array,k,L,p,minval));
    end
end



e1=plot(I_k_fromf,I_k_fromftilde,'ks');
e1.MarkerFaceColor='k';
e1.MarkerSize=8;
e1.LineWidth=2;
e1.HandleVisibility='off';
hold on;

C0=1;
f=@(x) C0*(x);
fun=fplot(f,[0. 1.],'LineWidth',2,'DisplayName','CFT','LineStyle','-.','LineWidth',2);
fun.Color='m';
set(get(get(fun,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on;

dim = [0.45 0.05 0.3 0.3];
str = {['$\varphi = \left(K/\zeta\right) = $' num2str(varphi)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',28,'LineStyle','-');


% dim = [0.3 0. 0.3 0.3];
% str = { 'NTRAJ : $5\times10^5$'};
% %str = {'$b=100$','$\epsilon=1.0,h^{*}=0.3$'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',22);



% h1=refline([0. r2_eq]);
% h1.LineWidth=3;
% h1.Color=[0. 0.5 0.];
% h1.DisplayName='Analytical';

% h2=line([20. 20.], [0. 10.],'Color','k','LineStyle','-.','LineWidth',3.);
% h2.HandleVisibility='off';

[h,icons,plots,legend_text]=legend({},'Location','northeast','FontSize',20,'Interpreter','latex','Box','off');

