% Total stress

% Create axes
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
% axes1.XScale='log';
axes1.YScale='log';
title('3D multi-bead-spring-dashpot','Interpreter','latex','FontSize',20);
% xlabel('$t^*$','FontSize',30,'Interpreter','latex');
xlabel('$1/N$','FontSize',30,'Interpreter','latex');
% y=ylabel('$\|{\mathbf{D}^{T}-\mathbf{D}}\|_2$','FontSize',42,'Interpreter','latex',...
%     'Rotation',90);
% set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlim([0. 0.2]);
% ylim([0. 1.]);
% pbaspect([1. 1. 1.]);
format long;
% % grid on;



%number of dimensions
ndim=3;






%ensemble_size
nsamp=200;
%\varphi=K/zeta
varphi=0;
dat=zeros(nsamp,5);

%%%%%%%%% set 1 %%%%%%%%%%%%%%%%%%%%

for i=1:nsamp
   
    %number of springs
    N=i;
    %size of diffusion matrix
    ndsize=N*ndim;
    [rb_mat] = construct_rouse_block_mat(N,ndim);

    %evaluating the diffusion matrix and its characteristics
    max_eig=max(eig(rb_mat));%storing the largest eigen-value
    min_eig=min(eig(rb_mat));%storing the smallest eigen-value
    devsym=rb_mat-rb_mat'; %deviation from symmetricity
    lin_devsym=reshape(devsym,[1 (ndsize*ndsize)]);
    normdev=norm(lin_devsym); %2-norm
    
    dat(i,1)=1./N;
    dat(i,2)=normdev;
    dat(i,3)=(1./min_eig);
    dat(i,4)=(max_eig/min_eig);
    dat(i,5)=cond(rb_mat);

end

xval=dat(:,1);
yval=dat(:,5);
xlog=log(xval);
ylog=log(yval);

% e1=plot(dat(:,1),dat(:,2),'ro');
% e1.MarkerFaceColor='r';
% e1.MarkerSize=10;
% e1.LineWidth=2;
% e1.DisplayName=['$\|{\mathbf{D}^{T}-\mathbf{D}}\|_2\,\varphi = $ ' num2str(varphi1)];
% hold on;
% 
% e2=plot(dat2(:,1),dat2(:,2),'ko');
% e2.MarkerFaceColor='k';
% e2.MarkerSize=10;
% e2.LineWidth=2;
% e2.DisplayName=['$\quad\quad\quad\quad\quad\varphi = $ ' num2str(varphi2)];
% hold on;
% 
% e3=plot(dat(:,1),dat(:,3),'rd');
% % e2.MarkerFaceColor='r';
% e3.MarkerSize=10;
% e3.LineWidth=2;
% e3.DisplayName=['$\lambda_{\mathrm{min}}^{-1}$'];
% hold on;
% 
% e4=plot(dat(:,1),dat(:,4),'bo');
% % e2.MarkerFaceColor='r';
% e4.MarkerSize=10;
% e4.LineWidth=2;
% e4.DisplayName=['$\lambda_{\mathrm{max}}/\lambda_{\mathrm{min}}$'];
% hold on;

% e5=plot(dat(:,1),dat(:,5),'ks');
% e5.MarkerSize=10;
% e5.LineWidth=2;
% e5.DisplayName=['$cond(\mathbf{D}^{T})$'];
% hold on;

e5=plot(xval,yval,'ks');
% e2.MarkerFaceColor='r';
e5.MarkerSize=10;
e5.LineWidth=2;
e5.DisplayName=['$cond(\mathbf{D}^{T})$'];
hold on;




% 
% e4=plot(dat2(:,1),dat2(:,3),'kd');
% % e2.MarkerFaceColor='r';
% e4.MarkerSize=10;
% e4.LineWidth=2;
% e4.DisplayName=['$\qquad\qquad\quad\varphi = $ ' num2str(varphi2)];
% hold on;

% C0=2.5;
% offset=-1.5;
% f=@(x) x.^2;
% % f=@(x) C0*(x)+offset;
% fun=fplot(f,[40. 150.],'LineWidth',2,'DisplayName','CFT','LineStyle','-.','LineWidth',2);
% fun.Color='r';
% set(get(get(fun,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% hold on;

dim = [0.45 0.1 0.3 0.3];
str = {['$\varphi = \left(K/\zeta\right) = $' num2str(varphi)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',24,'LineStyle','-');

% dim = [0.15 0.4 0.3 0.3];
% str = {'Red line has slope = 2'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',24,'LineStyle','-');


% dim = [0.45 0.1 0.3 0.3];
% str = {['$N = $' num2str(N)]};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',24,'LineStyle','-');




[h,icons,plots,legend_text]=legend({},'Location','east','FontSize',20,'Interpreter','latex','Box','on');

