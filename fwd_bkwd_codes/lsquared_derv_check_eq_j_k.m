
% Create axes
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
%axes1.XScale='log';
% axes1.YScale='log';
title('3D multi-bead-spring-dashpot','Interpreter','latex','FontSize',20);
% xlabel('$t^*$','FontSize',30,'Interpreter','latex');
xlabel('$Q_{j}^{(x)}$','FontSize',30,'Interpreter','latex');
% y=ylabel('$I_k$','FontSize',42,'Interpreter','latex',...
%     'Rotation',90);
% set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
% xlim([0. 5.]);
% ylim([0. 1.]);
pbaspect([1. 1. 1.]);
format long;
grid on;


%number of springs
N=10;
%number of dimensions
ndim=3;
%Value of "k". The "L_k" value that is evaluated
k=4;
%value of "j". Connector vector with respect to which gradient is measured
j=4;

rng(16032020);

%\omega=varphi/zeta
varphi=20;

p=(varphi/((2*varphi)+1))^2;

%creating initial configurations
Q=normrnd(0,1,[N,ndim]);
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%CALCULATING THE DERIVATIVE NUMERICALLY%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

del=1e-4;
nlim=80000;
yval=zeros(nlim,1);
yval(1)=L(k)*L(k);
xval=zeros(nlim,1);
% Q(j,1)=0.0;
xval(1)=Q(j,1);

for i=2:nlim
    Q(j,1) = Q(j,1)+del;
    xval(i)=Q(j,1);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    yval(i)=L(k)*L(k);
end

derv_complete=gradient(yval,del);

%calculation of derivative, point A
loc=10000;
xloc=xval(loc);
drvA=(yval(loc+1)-yval(loc-1))./(2*del);
configA=Q;
configA(j,1)=xloc;


%calculation of derivative, point B
loc=27000;
xloc=xval(loc);
drvB=(yval(loc+1)-yval(loc-1))./(2*del);
configB=Q;
configB(j,1)=xloc;

%calculation of derivative, point C
loc=60000;
xloc=xval(loc);
drvC=(yval(loc+1)-yval(loc-1))./(2*del);
configC=Q;
configC(j,1)=xloc;


e1=plot(xval,yval,'-b','DisplayName','$L^2_k$');
e1.MarkerFaceColor='b';
e1.MarkerSize=12;
e1.LineWidth=2;
% e1.HandleVisibility='off';
hold on;

e2=plot(xval,derv_complete,'-.r','DisplayName','$\partial L^2_k/\partial Q_{j}^{(x)}$');
e2.MarkerFaceColor='r';
e2.MarkerSize=12;
e2.LineWidth=2;
% e2.HandleVisibility='off';
hold on;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%CALCULATING THE DERIVATIVE ANALYTICALLY%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%% at point A %%%%%%%

Q=configA;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
y_at_A=L(k)*L(k);
drv_anltc_A = derv_lsq_k_q_k(L,Q,normQ,k,N);
chkA=2*L(k)*derv_l_k_q_k(L,Q,normQ,k,N);

p1=plot(Q(j,1),drv_anltc_A(1),'d','DisplayName','A');
p1.MarkerFaceColor=[0. 0.5 0.];
p1.Color=[0. 0.5 0.];
p1.MarkerSize=10;
p1.LineWidth=2;
p1.HandleVisibility='off';

% %%%%%% at point B %%%%%%%
% 
Q=configB;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
y_at_B=L(k)*L(k);
drv_anltc_B = derv_lsq_k_q_k(L,Q,normQ,k,N);

p2=plot(Q(j,1),drv_anltc_B(1),'o','DisplayName','B');
p2.MarkerFaceColor=[0. 0.5 0.];
p2.Color=[0. 0.5 0.];
p2.MarkerSize=10;
p2.LineWidth=2;
p2.HandleVisibility='off';
% 

%%%%%% at point C %%%%%%%

Q=configC;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
y_at_C=L(k)*L(k);
drv_anltc_C = derv_lsq_k_q_k(L,Q,normQ,k,N);

% 
p3=plot(Q(j,1),drv_anltc_C(1),'s','DisplayName','C');
p3.MarkerFaceColor=[0. 0.5 0.];
p3.Color=[0. 0.5 0.];
p3.MarkerSize=10;
p3.LineWidth=2;
p3.HandleVisibility='off';



dim = [0.55 0.05 0.3 0.3];
str = {['$\varphi = $' num2str(varphi)],['$\Delta = $' num2str(del)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'LineStyle','-');

dim = [0.55 0.40 0.3 0.3];
str = {['$N_{\mathrm{s}} = $' num2str(N)],['$k = $' num2str(k),',\,$j = $' num2str(j)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'LineStyle','-');



[h,icons,plots,legend_text]=legend({},'Location','northeast','FontSize',24,'Interpreter','latex','Box','on');

