% Create axes
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
%axes1.XScale='log';
%axes1.YScale='log';
title('3D multi-bead-spring-dashpot','Interpreter','latex','FontSize',20);
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
N=6;
%number of dimensions
ndim=3;
%Value of "i". The "M" value that is evaluated
i=2;
%Value of "k". The "P" value that is evaluated
k=2;
%value of "j". Connector vector with respect to which gradient is measured
j=3;

%Next we consider the term (D_m/D_n)
m=5;
n=3;

%Next we consider the tensor, T=(\bm{Q}_f\bm{Q}_s)/(Q_fQ_s)
f=3;
s=3;

% rng('shuffle');
rng(1687191421);

%\varphi=K/zeta
varphi=200;

p=(varphi/((2*varphi)+1))^2;

%creating initial configurations
Q=normrnd(0,1,[N,ndim]);
% Q(j,:)=0.05;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);

initQ=Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%CALCULATING THE DERIVATIVE NUMERICALLY%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
del=1e-6;
nlim=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FORWARD DIRECTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xval=zeros(nlim,1);
yval=zeros(nlim,1);
zval=zeros(nlim,1);

T11_fw=zeros(nlim,1);
T12_fw=zeros(nlim,1);
T13_fw=zeros(nlim,1);

T21_fw=zeros(nlim,1);
T22_fw=zeros(nlim,1);
T23_fw=zeros(nlim,1);

T31_fw=zeros(nlim,1);
T32_fw=zeros(nlim,1);
T33_fw=zeros(nlim,1);



% w.r.t x-direction

chdim=1;
xval(1)=Q(j,chdim);
temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
T11_fw(1)=temp(1,1);
T12_fw(1)=temp(1,2);
T13_fw(1)=temp(1,3);
for spcounter=2:nlim
    Q(j,chdim) = Q(j,chdim)+del;
    xval(spcounter)=Q(j,chdim);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
    T11_fw(spcounter)=temp(1,1);
    T12_fw(spcounter)=temp(1,2);
    T13_fw(spcounter)=temp(1,3);
end


%getting back to initial configuration
%before finding next set of gradients

Q=initQ;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);


% w.r.t y-direction

chdim=2;
yval(1)=Q(j,chdim);
temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
T21_fw(1)=temp(2,1);
T22_fw(1)=temp(2,2);
T23_fw(1)=temp(2,3);
for spcounter=2:nlim
    Q(j,chdim) = Q(j,chdim)+del;
    yval(spcounter)=Q(j,chdim);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
    T21_fw(spcounter)=temp(2,1);
    T22_fw(spcounter)=temp(2,2);
    T23_fw(spcounter)=temp(2,3);
end


%getting back to initial configuration
%before finding next set of gradients

Q=initQ;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);

% w.r.t z-direction

chdim=3;
zval(1)=Q(j,chdim);
temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
T31_fw(1)=temp(3,1);
T32_fw(1)=temp(3,2);
T33_fw(1)=temp(3,3);
for spcounter=2:nlim
    Q(j,chdim) = Q(j,chdim)+del;
    zval(spcounter)=Q(j,chdim);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
    T31_fw(spcounter)=temp(3,1);
    T32_fw(spcounter)=temp(3,2);
    T33_fw(spcounter)=temp(3,3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% BACKWARD DIRECTION %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q=initQ;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);

xval_bk=zeros(nlim,1);
yval_bk=zeros(nlim,1);
zval_bk=zeros(nlim,1);

T11_bk=zeros(nlim,1);
T12_bk=zeros(nlim,1);
T13_bk=zeros(nlim,1);

T21_bk=zeros(nlim,1);
T22_bk=zeros(nlim,1);
T23_bk=zeros(nlim,1);

T31_bk=zeros(nlim,1);
T32_bk=zeros(nlim,1);
T33_bk=zeros(nlim,1);



% w.r.t x-direction

chdim=1;
xval_bk(1)=Q(j,chdim);
temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
T11_bk(1)=temp(1,1);
T12_bk(1)=temp(1,2);
T13_bk(1)=temp(1,3);
for spcounter=2:nlim
    Q(j,chdim) = Q(j,chdim)-del;
    xval_bk(spcounter)=Q(j,chdim);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
    T11_bk(spcounter)=temp(1,1);
    T12_bk(spcounter)=temp(1,2);
    T13_bk(spcounter)=temp(1,3);
end


%getting back to initial configuration
%before finding next set of gradients

Q=initQ;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);


% w.r.t y-direction

chdim=2;
yval_bk(1)=Q(j,chdim);
temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
T21_bk(1)=temp(2,1);
T22_bk(1)=temp(2,2);
T23_bk(1)=temp(2,3);
for spcounter=2:nlim
    Q(j,chdim) = Q(j,chdim)-del;
    yval_bk(spcounter)=Q(j,chdim);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
    T21_bk(spcounter)=temp(2,1);
    T22_bk(spcounter)=temp(2,2);
    T23_bk(spcounter)=temp(2,3);
end


%getting back to initial configuration
%before finding next set of gradients

Q=initQ;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);

% w.r.t z-direction

chdim=3;
zval_bk(1)=Q(j,chdim);
temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
T31_bk(1)=temp(3,1);
T32_bk(1)=temp(3,2);
T33_bk(1)=temp(3,3);
for spcounter=2:nlim
    Q(j,chdim) = Q(j,chdim)-del;
    zval_bk(spcounter)=Q(j,chdim);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    temp=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
    T31_bk(spcounter)=temp(3,1);
    T32_bk(spcounter)=temp(3,2);
    T33_bk(spcounter)=temp(3,3);
end

T11_bk(1)=[];
T12_bk(1)=[];
T13_bk(1)=[];
T21_bk(1)=[];
T22_bk(1)=[];
T23_bk(1)=[];
T31_bk(1)=[];
T32_bk(1)=[];
T33_bk(1)=[];

T11=[flipud(T11_bk);T11_fw];
T12=[flipud(T12_bk);T12_fw];
T13=[flipud(T13_bk);T13_fw];


T21=[flipud(T21_bk);T21_fw];
T22=[flipud(T22_bk);T22_fw];
T23=[flipud(T23_bk);T23_fw];


T31=[flipud(T31_bk);T31_fw];
T32=[flipud(T32_bk);T32_fw];
T33=[flipud(T33_bk);T33_fw];

derv_T11=gradient(T11,del);
derv_T21=gradient(T21,del);
derv_T31=gradient(T31,del);

derv_T12=gradient(T12,del);
derv_T22=gradient(T22,del);
derv_T32=gradient(T32,del);

derv_T13=gradient(T13,del);
derv_T23=gradient(T23,del);
derv_T33=gradient(T33,del);
% 
derv_inbuilt_z=derv_T13+derv_T23+derv_T33;

%calculating x-component of divergence
divT_x_numeric=derv_T11+derv_T21+derv_T31;
%calculating y-component of divergence
divT_y_numeric=derv_T12+derv_T22+derv_T32;
%calculating z-component of divergence
divT_z_numeric=derv_T13+derv_T23+derv_T33;

divT_numeric=[divT_x_numeric divT_y_numeric divT_z_numeric];


% 
% e1=plot(xval,T11,'-','DisplayName','$f$');
% e1.Color=funcwheel(chdim);
% e1.MarkerFaceColor=funcwheel(chdim);
% e1.MarkerSize=12;
% e1.LineWidth=2;
% % e1.HandleVisibility='off';
% hold on;
% % 
% e2=plot(xval,derv_T11,'-.','DisplayName',['$ ' plug ' $']);
% e2.Color=funcwheel(chdim);
% e2.MarkerFaceColor=funcwheel(chdim);
% e2.MarkerSize=12;
% e2.LineWidth=2;
% % e2.HandleVisibility='off';
% hold on;

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%CALCULATING THE DERIVATIVE ANALYTICALLY%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % The formula is div(fT)=(grad(f)).T + f(div(T))
% % %
% % %

%%%%%% at point A %%%%%%%

Q=initQ;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
dyad=(Q(f,:)'*Q(s,:))./(normQ(f)*normQ(s));

M = fwd_coeff_all(p,L,N);
P = bkwd_coeff_all(p,L,N);

fac1=(1./(1.-M(i)-P(k)));
fac2=bkwd_coeff_poly(m,L,p,N)/bkwd_coeff_poly(n,L,p,N);
fac3=L(j)*L(j-1);

[derv_fac1,denom] = derv_pi_mk(i,k,j,p,L,Q,normQ,N);
derv_fac2=derv_ratio_d_m_n_q_j(m,n,j,N,L,p,Q,normQ);
derv_fac3=derv_l_kmin1_l_k_q_k(L,Q,normQ,j,N);

term1=fac1*fac2*derv_fac3;
term2=fac1*fac3*derv_fac2;
term3=fac2*fac3*derv_fac1;

%(grad(f)).T
rhs_t1=(term1+term2+term3)*dyad;

%f(div(T))
rhs_t2=(fac1*fac2*fac3)*derv_q_tens(f,s,j,Q,normQ);

drv_anltc_A=rhs_t1+rhs_t2;

residA=norm(divT_numeric(nlim,:)-drv_anltc_A)/norm(divT_numeric(nlim,:));

% 
% %from mixed formula
% [drv_anltc_A_mx,denom] = derv_pi_mk(i,k,j,p,L,Q,normQ,N);
% val_at_A_mx=1./denom;
% 
% %from backward substitution coefficient formula
% drv_anltc_A_bkwd = derv_ratio_d_m_n_q_j(k+1,k,j,N,L,p,Q,normQ);
% val_at_A_bkwd=bkwd_coeff_poly(k+1,L,p,N)/bkwd_coeff_poly(k,L,p,N);
% 
% residA=drv_anltc_A_mx-drv_anltc_A_bkwd;
% 
% p1=plot(Q(j,chdim),drv_anltc_A_mx(chdim),'d','DisplayName','A');
% p1.Color=pointwheel(chdim,:);
% p1.MarkerFaceColor=pointwheel(chdim,:);
% p1.MarkerSize=10;
% p1.LineWidth=2;
% p1.HandleVisibility='off';
% % % % 
% %%%%%%% at point B %%%%%%%
% 
% Q=configB;
% normQ = construct_norm(Q,N);
% L = constructL(Q,normQ,N);
% 
% %from mixed formula
% [drv_anltc_B_mx,denom] = derv_pi_mk(i,k,j,p,L,Q,normQ,N);
% val_at_B_mx=1./denom;
% 
% %from backward substitution coefficient formula
% drv_anltc_B_bkwd = derv_ratio_d_m_n_q_j(k+1,k,j,N,L,p,Q,normQ);
% val_at_B_bkwd=bkwd_coeff_poly(k+1,L,p,N)/bkwd_coeff_poly(k,L,p,N);
% 
% residB=drv_anltc_B_mx-drv_anltc_B_bkwd;
% 
% % 
% p2=plot(Q(j,chdim),drv_anltc_B_mx(chdim),'o','DisplayName','B');
% p2.Color=pointwheel(chdim,:);
% p2.MarkerFaceColor=pointwheel(chdim,:);
% p2.MarkerSize=10;
% p2.LineWidth=2;
% p2.HandleVisibility='off';
% % % 
% % %%%%%% at point C %%%%%%%
% % 
% Q=configC;
% normQ = construct_norm(Q,N);
% L = constructL(Q,normQ,N);
% 
% %from mixed formula
% [drv_anltc_C_mx,denom] = derv_pi_mk(i,k,j,p,L,Q,normQ,N);
% val_at_C_mx=1./denom;
% 
% %from backward substitution coefficient formula
% drv_anltc_C_bkwd = derv_ratio_d_m_n_q_j(k+1,k,j,N,L,p,Q,normQ);
% val_at_C_bkwd=bkwd_coeff_poly(k+1,L,p,N)/bkwd_coeff_poly(k,L,p,N);
% 
% residC=drv_anltc_C_mx-drv_anltc_C_bkwd;
% 
% 
% 
% p3=plot(Q(j,chdim),drv_anltc_C_mx(chdim),'s','DisplayName','C');
% p3.Color=pointwheel(chdim,:);
% p3.MarkerFaceColor=pointwheel(chdim,:);
% p3.MarkerSize=10;
% p3.LineWidth=2;
% p3.HandleVisibility='off';



% dim = [0.56 0.15 0.3 0.3];
% str = {['$\varphi = $' num2str(varphi)],['$\Delta = $' num2str(del)]};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'LineStyle','-');
% 
% 
% dim = [0.55 0.40 0.3 0.3];
% str = {['$N_{\mathrm{s}} = $' num2str(N)],['$i= $' num2str(i) ',\,$k = $' num2str(k),',\,$j = $' num2str(j)]};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'LineStyle','-');


xlabel(['$ ' xlabelstring ' $'],'FontSize',30,'Interpreter','latex');
[h,icons,plots,legend_text]=legend({},'Location','west','FontSize',20,'Interpreter','latex','Box','on');

