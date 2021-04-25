% Create axes
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
%axes1.XScale='log';
%axes1.YScale='log';
title('3D multi-bead-spring-dashpot, $f=1/(1-M_i-P_k)$','Interpreter','latex','FontSize',20);
% xlabel('$t^*$','FontSize',30,'Interpreter','latex');
% xlabel('$Q_{j}^{(x)}$','FontSize',30,'Interpreter','latex');
% y=ylabel('$1/\left(1-M_{k}\right)$','FontSize',42,'Interpreter','latex',...
%     'Rotation',90);
% set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlim([-1. 3.]);
% ylim([0. 1.]);
pbaspect([1. 1. 1.]);
format long;
grid on;

%dimension w.r.t which gradient is calculated
chdim=1;

cartesian=['x','y','z'];
funcwheel=['r','k','b'];
pointwheel=[0. 0.5 0.; 0.4940 0.1840 0.5560; 0.6350 0.0780 0.1840];

xs='Q_{j}^{(%s)}';
specstring='\\partial f/\\partial Q_{j}^{(%s)}';
xlabelstring=sprintf(xs,cartesian(chdim));
plug=sprintf(specstring,cartesian(chdim));

%number of springs
N=50;
%number of dimensions
ndim=3;
%Value of "i". The "M" value that is evaluated
i=12;
%Value of "k". The "P" value that is evaluated
k=12;
%value of "j". Connector vector with respect to which gradient is measured
j=11;

% rng(1652632735);
rng(1685304566);

%\varphi=K/zeta
varphi=200;

p=(varphi/((2*varphi)+1))^2;

%creating initial configurations
Q=normrnd(0,1,[N,ndim]);
% Q(j,chdim)=0.002;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
M = fwd_coeff_all(p,L,N);
M_i = fwd_coeff_specific(p,L,i);
P = bkwd_coeff_all(p,L,N);
P_k = bkwd_coeff_specific(p,L,k,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%CALCULATING THE DERIVATIVE NUMERICALLY%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

del=0.0001;
nlim=50000;
yval=zeros(nlim,1);
yval(1)=1./(1.-M_i-P_k);
xval=zeros(nlim,1);
xval(1)=Q(j,chdim);

for spcounter=2:nlim
    Q(j,chdim) = Q(j,chdim)+del;
    xval(spcounter)=Q(j,chdim);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    M_i = fwd_coeff_specific(p,L,i);
    P_k = bkwd_coeff_specific(p,L,k,N);
    yval(spcounter)=1./(1.-M_i-P_k);
end


%calculation of derivative, point A
loc=2050;
xloc=xval(loc);
drvA=(yval(loc+1)-yval(loc-1))./(2*del);

configA=Q;
configA(j,chdim)=xloc;

xloca=yval(loc);

%calculation of derivative, point B
loc=7300;
xloc=xval(loc);
drvB=(yval(loc+1)-yval(loc-1))./(2*del);


configB=Q;
configB(j,chdim)=xloc;

xlocb=yval(loc);

%calculation of derivative, point C
loc=30000;
xloc=xval(loc);
drvC=(yval(loc+1)-yval(loc-1))./(2*del);

configC=Q;
configC(j,chdim)=xloc;

xlocc=yval(loc);

derv_complete=gradient(yval,del);

e1=plot(xval,yval,'-','DisplayName','$f$');
e1.Color=funcwheel(chdim);
e1.MarkerFaceColor=funcwheel(chdim);
e1.MarkerSize=12;
e1.LineWidth=2;
% e1.HandleVisibility='off';
hold on;

e2=plot(xval,derv_complete,'-.','DisplayName',['$ ' plug ' $']);
e2.Color=funcwheel(chdim);
e2.MarkerFaceColor=funcwheel(chdim);
e2.MarkerSize=12;
e2.LineWidth=2;
% e2.HandleVisibility='off';
hold on;

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%CALCULATING THE DERIVATIVE ANALYTICALLY%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
%%%%%% at point A %%%%%%%

Q=configA;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);

[drv_anltc_A,denom] = derv_pi_mk(i,k,j,p,L,Q,normQ,N);
val_at_A=1./denom;


p1=plot(Q(j,chdim),drv_anltc_A(chdim),'d','DisplayName','A');
p1.Color=pointwheel(chdim,:);
p1.MarkerFaceColor=pointwheel(chdim,:);
p1.MarkerSize=10;
p1.LineWidth=2;
p1.HandleVisibility='off';
% % % % 
%%%%%%% at point B %%%%%%%

Q=configB;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);

[drv_anltc_B,denom] = derv_pi_mk(i,k,j,p,L,Q,normQ,N);
val_at_B=1./denom;

% 
p2=plot(Q(j,chdim),drv_anltc_B(chdim),'o','DisplayName','B');
p2.Color=pointwheel(chdim,:);
p2.MarkerFaceColor=pointwheel(chdim,:);
p2.MarkerSize=10;
p2.LineWidth=2;
p2.HandleVisibility='off';
% % 
% %%%%%% at point C %%%%%%%
% 
Q=configC;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
[drv_anltc_C,denom] = derv_pi_mk(i,k,j,p,L,Q,normQ,N);
val_at_C=1./denom;

p3=plot(Q(j,chdim),drv_anltc_C(chdim),'s','DisplayName','C');
p3.Color=pointwheel(chdim,:);
p3.MarkerFaceColor=pointwheel(chdim,:);
p3.MarkerSize=10;
p3.LineWidth=2;
p3.HandleVisibility='off';



dim = [0.56 0.05 0.3 0.3];
str = {['$\varphi = $' num2str(varphi)],['$\Delta = $' num2str(del)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'LineStyle','-');


dim = [0.55 0.40 0.3 0.3];
str = {['$N_{\mathrm{s}} = $' num2str(N)],['$i= $' num2str(i) ',\,$k = $' num2str(k),',\,$j = $' num2str(j)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'LineStyle','-');


xlabel(['$ ' xlabelstring ' $'],'FontSize',30,'Interpreter','latex');
[h,icons,plots,legend_text]=legend({},'Location','west','FontSize',20,'Interpreter','latex','Box','on');

