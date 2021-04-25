% Total stress

% Create axes
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
%axes1.XScale='log';
%axes1.YScale='log';
title('3D multi-bead-spring-dashpot, $N=30$','Interpreter','latex','FontSize',20);
xlabel('$k$','FontSize',30,'Interpreter','latex');
y=ylabel('$I_k$','FontSize',42,'Interpreter','latex',...
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
% %Value of "k". The "M" value that is evaluated
% k=8;
% %value of "j". Connector vector with respect to which gradient is measured
% j=7;

rng(16032020);

%\omega=K/zeta
omega=200;

p=(omega/((2*omega)+1))^2;

%creating initial configurations
Q=normrnd(0,1,[N,ndim]);
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
M = fwd_coeff_all(p,L,N);

k = [1:N];
I_k_fromf = zeros(N,1);

for i=1:N
    I_k_fromf(i)=fwd_coeff_poly(i,L,p,N);
end



e1=plot(k,I_k_fromf,'ro');
e1.MarkerFaceColor='r';
e1.MarkerSize=8;
e1.LineWidth=2;
e1.HandleVisibility='off';
hold on;


dim = [0.52 0.6 0.3 0.3];
str = {['$\omega = \left(K/\zeta\right) = $' num2str(omega)],['$p = $' num2str(p)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',22,'LineStyle','-');

dim = [0.52 0.4 0.3 0.3];
str = {'$I_{k}=\prod_{i=1}^{k}\left(1-M_i\right)$',...
    '$M_{k}=p\left(\frac{L^2_{k-1}}{1-M_{k-1}}\right)$','$M_1=0$'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',22,'LineStyle','-');


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

