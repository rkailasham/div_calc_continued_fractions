function [drvA,drvB,drvC,drv_anltc_A,drv_anltc_B,drv_anltc_C] = d_k_derv_anltc_numeric(chdim,N,Q,j,k,L,varphi,del,nlim)

cartesian=['x','y','z'];
funcwheel=['r','k','b'];
pointwheel=[0. 0.5 0.; 0.4940 0.1840 0.5560; 0.6350 0.0780 0.1840];

xs='Q_{j}^{(%s)}';
specstring='\\partial D_k/\\partial Q_{j}^{(%s)}';
xlabelstring=sprintf(xs,cartesian(chdim));
plug=sprintf(specstring,cartesian(chdim));

p=(varphi/((2*varphi)+1))^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%CALCULATING THE DERIVATIVE NUMERICALLY%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


yval=zeros(nlim,1);
xval=zeros(nlim,1);

xval(1)=Q(j,chdim);
yval(1)=bkwd_coeff_poly(k,L,p,N);



for i=2:nlim
    Q(j,chdim) = Q(j,chdim)+del;
    xval(i)=Q(j,chdim);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    yval(i)=bkwd_coeff_poly(k,L,p,N);
end


derv_complete=gradient(yval,del);
%calculation of derivative, point A
loc=4600;
xloc=xval(loc);
drvA=(yval(loc+1)-yval(loc-1))./(2*del);
configA=Q;
configA(j,chdim)=xloc;


%calculation of derivative, point B
loc=21500;
xloc=xval(loc);
drvB=(yval(loc+1)-yval(loc-1))./(2*del);
configB=Q;
configB(j,chdim)=xloc;

%calculation of derivative, point C
loc=38000;
xloc=xval(loc);
drvC=(yval(loc+1)-yval(loc-1))./(2*del);
configC=Q;
configC(j,chdim)=xloc;
% % 


e1=plot(xval,yval,'-','DisplayName','$D_k$');
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
% %%%%%% at point A %%%%%%%
% 
Q=configA;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
D_k_at_A=fwd_coeff_poly(k,L,p,N);
[drv_anltc_A] = derv_d_k_q_j(k,j,L,p,Q,normQ,N);

p1=plot(Q(j,chdim),drv_anltc_A(chdim),'d','DisplayName','A');
p1.Color=pointwheel(chdim,:);
p1.MarkerFaceColor=pointwheel(chdim,:);
p1.MarkerSize=10;
p1.LineWidth=2;
p1.HandleVisibility='off';
% % 
%%%%%%% at point B %%%%%%%

Q=configB;
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
D_k_at_B=fwd_coeff_poly(k,L,p,N);
[drv_anltc_B] = derv_d_k_q_j(k,j,L,p,Q,normQ,N);
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
D_k_at_C=fwd_coeff_poly(k,L,p,N);
[drv_anltc_C] = derv_d_k_q_j(k,j,L,p,Q,normQ,N);
% 
p3=plot(Q(j,chdim),drv_anltc_C(chdim),'s','DisplayName','C');
p3.Color=pointwheel(chdim,:);
p3.MarkerFaceColor=pointwheel(chdim,:);
p3.MarkerSize=10;
p3.LineWidth=2;
p3.HandleVisibility='off';


xlabel(['$ ' xlabelstring ' $'],'FontSize',30,'Interpreter','latex');
[h,icons,plots,legend_text]=legend({},'Location','northeast','FontSize',20,'Interpreter','latex','Box','on');

end