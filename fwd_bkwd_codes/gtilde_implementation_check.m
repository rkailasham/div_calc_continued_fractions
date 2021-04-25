% Total stress

% Create axes
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
%axes1.XScale='log';
%axes1.YScale='log';
title('$D_k$ for 3D multi-bead-spring-dashpot, $N=30$','Interpreter','latex','FontSize',20);
% xlabel('$t^*$','FontSize',30,'Interpreter','latex');
xlabel('Calculated from $g$','FontSize',30,'Interpreter','latex');
y=ylabel('Calculated from $\tilde{g}$','FontSize',36,'Interpreter','latex',...
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
%Value of "k". The "P" value that is evaluated
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
%%%%%  CALCULATING D_k using "g"            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D_k_fromg = zeros(N,1);

for k=1:N
    D_k_fromg(k)=bkwd_coeff_poly(k,L,p,N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  CALCULATING D_k using "gtilde"       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_k_fromgtilde = ones(N,1);

for k=1:N
    for nu=k:(N-1)
        compare_array=gen_list(nu); %this involves repetitive calculation
                                    %room for speeding up
        minval=min(nu);
        D_k_fromgtilde(k)=D_k_fromgtilde(k)+(L(nu)*L(nu)*bkwd_coeff_lsquared_nu_composite(nu,compare_array,k,L,p,minval,N));
    end
end

% 
% 
e1=plot(D_k_fromg,D_k_fromgtilde,'bd');
e1.MarkerFaceColor='b';
e1.MarkerSize=8;
e1.LineWidth=2;
e1.HandleVisibility='off';
hold on;

C0=1;
f=@(x) C0*(x);
fun=fplot(f,[0. 1.],'LineWidth',2,'DisplayName','CFT','LineStyle','-.','LineWidth',2);
fun.Color='r';
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

