
%number of springs
N=30;
%number of dimensions
ndim=3;
%reference index
k=5;


rng(40052020);

%\varphi=K/zeta
varphi=200;

p=(varphi/((2*varphi)+1))^2;

%creating initial configurations
Q=normrnd(0,1,[N,ndim]);
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
M = fwd_coeff_all(p,L,N);
P = bkwd_coeff_all(p,L,N);
dir_calc=1./(1.-M-P);


ref_dir=dir_calc(k);
ratio_calc=calc_inv_pi_mk(k,k,L,p,N);
diff=(ref_dir-ratio_calc);
perc_diff=abs(diff/ref_dir)*100.;

fin_dat=[N,k,ref_dir,ratio_calc,diff,perc_diff];
dlmwrite('output.dat', fin_dat,'precision','%20.15f','delimiter',' ');
