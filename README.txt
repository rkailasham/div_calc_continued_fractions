25-APR-2021: 

This repository contains the MATLAB scripts for the analytical calculation of divergences 
of finite continued fractions. The relevant formulae have been provided in the 
"Supplementary Material" section of the paper titled "Rouse model with fluctuating
internal friction". 


Given below is a list of scripts and their purpose. The equation numbers are 
specified with respect to the Supplementary Material section, unless
mentioned otherwise.


Within the folder "fwd_bkwd_codes",
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A: Forward continued product			     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fwd_coeff_poly.m			: Constructs I_k using Eq. (19)
f_rec.m					: Implements recursion given by Eq. (20)
fwd_coeff_lsquared_nu_composite.m	: Returns coefficients of L^2_{\nu} in Eq. (22)
qnu_k_tilde_coeff.m			: For Eq. (23)
ftilde_rec.m				: Implements recursion for "f_tilde" needed in Eq. (23)
ftilde_implementation_check.m		: Checks that I_k calculated from Eq. (19) and
					  Eq. (22) are identical
derv_qnu_k_tilde_q_j.m			: Implementation of Eq. (26)
derv_i_k_q_j.m				: Eq. (29)
driver_fwd_i_k_derv_check.m		: Generates Fig. S1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B: Backward continued product			     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
bkwd_coeff_poly.m			: Constructs D_k using Eq. (30)
g_rec.m					: Implements recursion given by Eq. (31)
bkwd_coeff_lsquared_nu_composite.m	: Returns coefficients of L^2_{\nu} in Eq. (32)
bkwd_qnu_k_tilde_coeff.m		: For Eq. (33)
gtilde_rec.m				: Implements recursion for "g_tilde" needed in Eq. (34)
gtilde_implementation_check.m		: Checks that D_k calculated from Eq. (30) and
                                          Eq. (32) are identical
bkwd_qnu_k_tilde_coeff.m		: Implementation of Eq. (36)
derv_d_k_q_j.m				: Eq. (37)
driver_bkwd_d_k_derv_check.m		: Generates Fig. S2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C: List of tensor identities			     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

derv_q_tens.m				: Eq. (40)
derv_l_kmin1_q_k.m			: Eq. (42)
derv_l_k_q_k.m				: Eq. (43)
derv_lsq_kmin1_q_k.m			: Eq. (45)
derv_lsq_k_q_k.m			: Eq. (46)
derv_l_kmin1_l_k_q_k.m			: Eq. (48)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F: Calculation of square root of the diffusion tensor			     			     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver_diffmat_sym_eig_direct.m         : Generates Fig. S6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E: Calculation of divergence terms in SDE and stress tensor expression			     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The folder "fig_s3_derv_calc" contains the necessary data for generating Fig. S3.
Within this folder,

driver_plot_derv_accuracy.m		: Driver routine that generates Fig. S3
derv_pi_mk.m				: Code that analytically evaluates the gradient in Fig. S3
					  This code was run on the various test cases, in order to 
					  generate data that is used for plotting Fig. S3
	
The folder "divergence_V_jk_trans" contains the data needed for plotting Fig. S4(a) and Fig. S5(a).
Within this folder,

div_accuracy.m				: Driver routine that generates Fig. S4(a)
div_timing.m				: Driver routine that generates Fig. S5(a)
divergence_quantify_error		: Code that analytically evaluates the divergence of $\bm{V}^{T}_{jk}$
                                          This code was run on the various test cases, in order to
                                          generate data used for plotting Fig. S4(b) and Fig. S5(b)

The folder "divergence_mu_kl_trans" contains the data needed for plotting Fig. S4(b) and Fig. S5(b).
Within this folder,

mu_div_accuracy.m			: Driver routine that generates Fig. S4(b)
mu_div_timing.m				: Driver routine that generates Fig. S5(b)
mu_divergence_quantify_error.m		: Code that analytically evaluates the divergence of $\bm{\mu}^{T}_{kl}$
					  This code was run on the various test cases, in order to 
					  generate data used for plotting Fig. S4(b) and Fig. S5(b)

