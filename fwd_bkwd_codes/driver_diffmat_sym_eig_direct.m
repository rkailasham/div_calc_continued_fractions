% Total stress

% Create axes
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
% axes1.XScale='log';
axes1.YScale='log';
% title('3D multi-bead-spring-dashpot','Interpreter','latex','FontSize',20);
% xlabel('$t^*$','FontSize',30,'Interpreter','latex');
xlabel('Random Sample No.','FontSize',30,'Interpreter','latex');
% y=ylabel('$\|{\mathbf{D}^{T}-\mathbf{D}}\|_2$','FontSize',42,'Interpreter','latex',...
%     'Rotation',90);
% set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
% xlim([0. 5.]);
ylim([1e-16 1.]);
% pbaspect([1. 1. 1.]);
format long;
% % grid on;


%number of springs
N=45;
%number of dimensions
ndim=3;
%size of diffusion matrix
ndsize=N*ndim;




% rng(27042021);


%ensemble_size
nsamp=100;

dat1=zeros(nsamp,3);
dat2=zeros(nsamp,3);
dat3=zeros(nsamp,3);

%checking for imaginary eigen values
imag_flag=0;
imag_count=0;
imag_count2=0;

%checking for negative eigen values
neg_flag=0;
neg_count=0;
neg_count2=0;

%%%%%%%%% set 1 %%%%%%%%%%%%%%%%%%%%
%\varphi=K/zeta
varphi1=5;
p=(varphi1/((2*varphi1)+1))^2;

for i=1:nsamp
    imag_flag=0;
    neg_flag=0;
    %creating initial configurations, shuffling the seed each time
    rng('shuffle');

    Q=normrnd(0,1,[N,ndim]);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);

    %evaluating the diffusion matrix and its characteristics
    [diffMat] = diffMat_eval_direct(varphi1,L,Q,normQ,N,ndim);
    diffMat_t=diffMat';
    min_eig=min(eig(diffMat_t));%storing the smallest eigen-value
    devsym=diffMat_t-diffMat_t'; %deviation from symmetricity
    lin_devsym=reshape(devsym,[1 (ndsize*ndsize)]);
    normdev=norm(lin_devsym); %2-norm
    list_eig=eig(diffMat_t);
    imag_list=imag(list_eig);
    chk_list=imag_list(imag_list~=0);
    neg_list=list_eig(list_eig<0);

    if(isempty(chk_list))
        imag_flag=0;
    else
        imag_flag=1;
        imag_count=imag_count+1;
    end
    
    if(isempty(neg_list))
        neg_flag=0;
    else
        neg_flag=1;
        neg_count=neg_count+1;
    end
    
    
    dat1(i,1)=i;
    dat1(i,2)=normdev;
    dat1(i,3)=min_eig;

end


%%%%%%%%% set 2 %%%%%%%%%%%%%%%%%%%%
%\varphi=K/zeta
varphi2=500;
p=(varphi2/((2*varphi2)+1))^2;

for i=1:nsamp
    imag_flag=0;
    %creating initial configurations, shuffling the seed each time
    rng('shuffle');

    Q=normrnd(0,1,[N,ndim]);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);

    %evaluating the diffusion matrix and its characteristics
    [diffMat] = diffMat_eval_direct(varphi2,L,Q,normQ,N,ndim);
    diffMat_t=diffMat';
    min_eig=min(eig(diffMat_t));%storing the smallest eigen-value
    devsym=diffMat_t-diffMat_t'; %deviation from symmetricity
    lin_devsym=reshape(devsym,[1 (ndsize*ndsize)]);
    normdev=norm(lin_devsym); %2-norm
    chk_list=imag_list(imag_list~=0);
    neg_list=list_eig(list_eig<0);

    if(isempty(chk_list))
        imag_flag=0;
    else
        imag_flag=1;
        imag_count2=imag_count2+1;
    end        
    
    if(isempty(neg_list))
        neg_flag=0;
    else
        neg_flag=1;
        neg_count2=neg_count2+1;
    end
    
    dat2(i,1)=i;
    dat2(i,2)=normdev;
    dat2(i,3)=min_eig;

end


e1=plot(dat1(:,1),dat1(:,2),'ro');
e1.MarkerFaceColor='r';
e1.MarkerSize=10;
e1.LineWidth=2;
e1.DisplayName=['{$|$\boldmath$\widehat{d}$}$|,\,\varphi = $ ' num2str(varphi1)];
hold on;

e2=plot(dat2(:,1),dat2(:,2),'ko');
e2.MarkerFaceColor='k';
e2.MarkerSize=10;
e2.LineWidth=2;
e2.DisplayName=['{$|$\boldmath$\widehat{d}$}$|,\,\varphi = $ ' num2str(varphi2)];
% e2.DisplayName=['$\|{\mathbf{D}^{T}-\mathbf{D}}\|\,\varphi = $ ' num2str(varphi2)];
hold on;

e3=plot(dat1(:,1),dat1(:,3),'rd');
% e2.MarkerFaceColor='r';
e3.MarkerSize=10;
e3.LineWidth=2;
e3.DisplayName=['$\lambda_{\mathrm{min}},\,\varphi = $ ' num2str(varphi1)];
% e3.DisplayName=['$\lambda_{\mathrm{min}}\quad\quad\quad\,\varphi = $ ' num2str(varphi1)];
hold on;

e4=plot(dat2(:,1),dat2(:,3),'kd');
% e2.MarkerFaceColor='r';
e4.MarkerSize=10;
e4.LineWidth=2;
e4.DisplayName=['$\lambda_{\mathrm{min}},\,\varphi = $ ' num2str(varphi2)];
% e4.DisplayName=['$\qquad\qquad\quad\varphi = $ ' num2str(varphi2)];
hold on;

% C0=1;
% f=@(x) C0*(x);
% fun=fplot(f,[1.0 1.5],'LineWidth',2,'DisplayName','CFT','LineStyle','-.','LineWidth',2);
% fun.Color='b';
% set(get(get(fun,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% hold on;

% dim = [0.45 0.1 0.3 0.3];
% str = {['$\varphi = \left(K/\zeta\right) = $' num2str(varphi)]};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',24,'LineStyle','-');


dim = [0.15 0.25 0.3 0.3];
str = {['$N = $' num2str(N)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',24,'LineStyle','-');




[h,icons,plots,legend_text]=legend({},'Location','east','FontSize',20,'Interpreter','latex','Box','on');

