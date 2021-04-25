% Create axes
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
%axes1.XScale='log';
% axes1.YScale='log';
% title('3D multi-bead-spring-dashpot','Interpreter','latex','FontSize',20);
% xlabel('$t^*$','FontSize',30,'Interpreter','latex');
% xlabel('$Q_{j}^{(x)}$','FontSize',30,'Interpreter','latex');
% y=ylabel('$I_k$','FontSize',42,'Interpreter','latex',...
%     'Rotation',90);
% set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlim([-1. 4.5]);
% ylim([0. 1.]);
% pbaspect([1. 1. 1.]);
format long;
grid on;



%number of springs
N=20;
%number of dimensions
ndim=3;
%Value of "k". The "I_k" value that is evaluated
k=14;
%value of "j". Connector vector with respect to which gradient is measured
j=7;

rng(722958551);
% rng('shuffle');

drvA=zeros(3,1);
drvB=zeros(3,1);
drvC=zeros(3,1);

drv_anltc_A=zeros(3,3);
drv_anltc_B=zeros(3,3);
drv_anltc_C=zeros(3,3);

%\varphi=K/zeta
varphi=20;

p=(varphi/((2*varphi)+1))^2;

%creating initial configurations
Q=normrnd(0,1,[N,ndim]);
% Q(j-1,:)=0.2;
% Q(j,:)=0.5;
% Q(j+1,:)=0.3;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);

del=0.0001;
nlim=50000;

%choose dimension to vary
% chdim=1;
% [drvA,drvB,drvC,drv_anltc_A,drv_anltc_B,drv_anltc_C]=i_k_derv_anltc_numeric(chdim,N,Q,j,k,L,varphi,del,nlim);

% [drvA(chdim),drvB(chdim),drvC(chdim),drv_anltc_A(chdim,:),drv_anltc_B(chdim,:),drv_anltc_C(chdim,:)]=i_k_derv_anltc_numeric(chdim,N,Q,j,k,L,varphi,del,nlim);
% chdim=2;
% [drvA(chdim),drvB(chdim),drvC(chdim),drv_anltc_A(chdim,:),drv_anltc_B(chdim,:),drv_anltc_C(chdim,:)]=i_k_derv_anltc_numeric(chdim,N,Q,j,k,L,varphi,del,nlim);
chdim=2;
[drvA(chdim),drvB(chdim),drvC(chdim),drv_anltc_A(chdim,:),drv_anltc_B(chdim,:),drv_anltc_C(chdim,:)]=i_k_derv_anltc_numeric(chdim,N,Q,j,k,L,varphi,del,nlim);

dim = [0.55 0.25 0.3 0.3];
str = {['$\varphi = $' num2str(varphi)],['$\Delta_{\mathrm{d}} = $' num2str(del)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'LineStyle','-');

dim = [0.55 0.40 0.3 0.3];
str = {['$N = $' num2str(N)],['$k = $' num2str(k),',\,$j = $' num2str(j)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'LineStyle','-');


% [h,icons,plots,legend_text]=legend({},'Location','northwest','FontSize',20,'Interpreter','latex','Box','on');

