
gen_data = function(n = 100, d = 2L, num_grids = 64L, grids_lim = c(-1,1), random = FALSE, poly_degree = 8L,
							a = 0.1, b = 20, center = NULL, rate = NULL, max_range = 6, thres = thres){
	## return a list containing all the variables we used to generate Y 
	## n: the number of subjects
	## d: dimension of the figure
	## num_grids: the grids we have on each dimension
	set.seed(2017)
 
	grids = GP.generate.grids(d = d, num_grids = num_grids, grids_lim = grids_lim)
	Xmat= GP.eigen.funcs.fast(grids = grids, poly_degree=poly_degree, a=a, b=b)
	lambda = GP.eigen.value(poly_degree=poly_degree, a=a, b=b, d=d)

	tau_xi_sq = 1
	c_l = rnorm(length(lambda), mean=0, sd=sqrt(lambda))
	xi = sqrt(tau_xi_sq) * (Xmat%*%c_l)
	## generate sigma1 and sigma2

	sigma_pos_sq = ifelse(xi > thres, xi, 0)
  sigma_neg_sq = ifelse(xi < -thres, -xi, 0)
	sigma_pos = sqrt(sigma_pos_sq)
  sigma_neg = sqrt(sigma_neg_sq)
	
	log_tau_1_sq = Xmat%*%rnorm(length(lambda), mean=0, sd=0.5*sqrt(lambda))

    log_tau_2_sq = Xmat%*%rnorm(length(lambda), mean=0, sd=0.5*sqrt(lambda))
   	tau_1_sq = 0.25*c(exp(log_tau_1_sq))
   	tau_2_sq = 0.25*c(exp(log_tau_2_sq))
	rho = (sigma_pos^2 - sigma_neg^2)/(sqrt(sigma_pos^2 + sigma_neg^2 + tau_1_sq) * sqrt(sigma_pos^2 + sigma_neg^2 + tau_2_sq))
	
	eta_pos = eta_neg = matrix(nrow = num_grids^d, ncol = n)
	e_pos = matrix(nrow = length(lambda), ncol = n)
	e_neg = matrix(nrow = length(lambda), ncol = n)
	for(i in c(1:n)){
		 e_pos[,i] = rnorm(length(lambda),mean=0,sd=sqrt(lambda))
		 e_neg[,i] = rnorm(length(lambda),mean=0,sd=sqrt(lambda))
	}
	
	eta_pos = c(sigma_pos) * (Xmat%*%e_pos)
	eta_neg = c(sigma_neg) * (Xmat%*%e_neg)
	
	eps_1 = eps_2 = matrix(nrow = num_grids^d, ncol = n)
	for(i in c(1:n)){
		eps_1[,i] = rnorm(length(tau_1_sq), 0, sd = sqrt(tau_1_sq))
    	eps_2[,i] = rnorm(length(tau_2_sq), 0, sd = sqrt(tau_2_sq))
	}
	Y_1 = Y_2 = matrix(nrow = num_grids^d, ncol = n)
	for(i in c(1:n)){
		Y_1[,i] = eta_pos[,i] + eta_neg[,i] + eps_1[,i]
		Y_2[,i] = eta_pos[,i] - eta_neg[,i] + eps_2[,i]
	}
	thres_xi = sigma_pos_sq - sigma_neg_sq
  	abs_thres_xi = abs(thres_xi)
  	var_Y_1 = tau_1_sq + abs_thres_xi
  	var_Y_2 = tau_2_sq + abs_thres_xi
	ratio = (sigma_pos^2 + sigma_neg^2)/(sigma_pos^2 + sigma_neg^2 + tau_1_sq^2 + tau_2_sq^2)
	return(list('x' = grids, 'true_c' = c_l, 'true_e_pos' = e_pos, 'true_e_neg' = e_neg, 'xi' = xi, 'Xmat' = Xmat,'sigma_pos_sq' = sigma_pos^2, 'sigma_neg_sq' = sigma_neg^2, 'tau_1_sq' = tau_1_sq,
	 'tau_2_sq' = tau_2_sq, 'tau_xi_sq' = tau_xi_sq, 'rho' = rho, 'eta_pos' = eta_pos, 'eta_neg' = eta_neg, 'eps_1' = eps_1, 'eps_2' = eps_2,
	 'Y_1' = Y_1, 'Y_2' = Y_2, 'var_Y_1' = var_Y_1, 'var_Y_2' = var_Y_2, 'ratio' = ratio))
}



gen_data_design = function(n = 50, d = 2L, num_grids = 64L, grids_lim = c(0,1), random = FALSE, poly_degree = 20L,
							a=0.1, b=20, tau_scale = 0.5, pos_radius = 0.2, neg_radius = 0.2, pos_mag = 0.75, neg_mag = 0.85){
    
    set.seed(2017)


	  grids = GP.generate.grids(d = d, num_grids = num_grids, grids_lim = grids_lim)
    Xmat= GP.eigen.funcs.fast(grids = grids, poly_degree=poly_degree, a=a, b=b)
	  lambda = GP.eigen.value(poly_degree=poly_degree, a=a, b=b, d=2)
  	log_tau_1_sq = Xmat%*%rnorm(length(lambda),mean=0,sd=tau_scale*sqrt(lambda))
  	log_tau_2_sq = Xmat%*%rnorm(length(lambda),mean=0,sd=tau_scale*sqrt(lambda))
  	tau_1_sq = 0.25*exp(log_tau_1_sq)
  	tau_2_sq = 0.25*exp(log_tau_2_sq)
  	#tau_1_sq = rinvgamma(nrow(grids), 1, 0.1)
  	#tau_2_sq = rinvgamma(nrow(grids), 1, 0.1)
    ############### package data
    # sigma_pos_sq = ifelse(abs(grids[,1]-0.3) + abs(grids[,2]-0.7) < 0.2, pos_mag, 0)
    # sigma_neg_sq = ifelse(abs(grids[,1]-0.7) + abs(grids[,2]-0.3) < 0.2, neg_mag, 0)

  	sigma_pos_sq = ifelse(abs(grids[,1]-0.3) + abs(grids[,2]-0.7) < pos_radius, pos_mag, 0)
  	sigma_pos_sq = ifelse(abs(grids[,1]-0.3) + abs(grids[,2]-0.3) < pos_radius, pos_mag, sigma_pos_sq)
  	sigma_neg_sq = ifelse((grids[,1]-0.5)^2 + (grids[,2]-0.5)^2 < neg_radius^2, neg_mag, 0)
    sigma_neg_sq = ifelse(abs(grids[,1]-0.7) + abs(grids[,2]-0.3) < neg_radius, neg_mag, sigma_neg_sq)
    sigma_pos_sq = ifelse(abs(grids[,1]-0.7) + abs(grids[,2]-0.7) < pos_radius, pos_mag, sigma_pos_sq)

  	sigma_pos = sqrt(sigma_pos_sq)
  	sigma_neg = sqrt(sigma_neg_sq)
  	eta_pos = matrix(NA, nrow=nrow(grids), ncol = n)
  	eta_neg = matrix(NA, nrow=nrow(grids), ncol = n)
  	Y_1 = matrix(NA, nrow=nrow(grids), ncol = n)
  	Y_2 = matrix(NA, nrow=nrow(grids), ncol = n)
  	for(i in 1:n){
    	eta_pos[,i] = as.numeric(sigma_pos) * (Xmat %*% rnorm(length(lambda),mean = 0,sd = sqrt(lambda)))
    	eta_neg[,i] = as.numeric(sigma_neg) * (Xmat %*% rnorm(length(lambda),mean = 0,sd = sqrt(lambda)))
    	
  	}
  	#set.seed(NULL)
  	
  	#tau_1_sq = rep(1, nrow(grids))
  	#tau_2_sq = rep(1, nrow(grids))
  	for(i in 1:n){
  		Y_1[,i] = rnorm(length(tau_1_sq), sd = sqrt(tau_1_sq))
    	Y_2[,i] = rnorm(length(tau_2_sq), sd = sqrt(tau_2_sq))
  	}
  	Y_1 = Y_1 + eta_pos + eta_neg
  	Y_2 = Y_2 + eta_pos - eta_neg
  	thres_xi = sigma_pos_sq - sigma_neg_sq
  	abs_thres_xi = abs(thres_xi)
  	var_Y_1 = tau_1_sq + abs_thres_xi
  	var_Y_2 = tau_2_sq + abs_thres_xi
  	rho = thres_xi/sqrt(var_Y_1 * var_Y_2)
  
  	cor_type = ifelse(rho > 0, 1, 0)
  	cor_type = ifelse(rho < 0, -1, cor_type)
  

	ratio = (sigma_pos^2 + sigma_neg^2)/(sigma_pos^2 + sigma_neg^2 + tau_1_sq^2 + tau_2_sq^2) #### signal to noise ratio
	return(list('x' = grids, 'Xmat' = Xmat,'sigma_pos_sq' = sigma_pos_sq, 'sigma_neg_sq' = sigma_neg_sq, 'tau_1_sq' = tau_1_sq,
	 	'tau_2_sq' = tau_2_sq, 'rho' = rho, 'eta_pos' = eta_pos, 'eta_neg' = eta_neg,
	   'Y_1' = Y_1, 'Y_2' = Y_2, 'var_Y_1' = var_Y_1, 'var_Y_2' = var_Y_2, 'ratio' = ratio, 'cor_type' = cor_type))
}
	



gen_data_3d = function(region_num, T, burn_in, d = 3L){
    ## My mac
    load(file = paste('/Users/moyan/Dropbox\ (University\ of\ Michigan)/Moyan/HCP/res_region_',region_num,'.Rdata', sep = ''))
    
    ## jinming
    #load(file = paste('/Users/neversalar/Moyan/Dropbox\ (University\ of\ Michigan)/Moyan/HCP/res_region_',region_num,'.Rdata', sep = ''))
    chain = res$chain
    Xmat = res$Xmat
    lambda = res$lambda
    filename = paste('/Users/moyan/Dropbox\ (University\ of\ Michigan)/Moyan-Jian-Share/data/residual.region',region_num,'.Rdata', sep = '')
    #filename = paste('/Users/neversalar/Moyan/Dropbox\ (University\ of\ Michigan)/Moyan-Jian-Share/data/residual.region',region_num,'.Rdata', sep = '')

    load(file = filename)
    grids = dat$grids
    index = c(as.integer(rownames(dat$grids)))
    #### added for new grids generation function
    # grids_list = list()
    # grids_list[[1]] = c(1:91)
    # grids_list[[2]] = c(1:109)
    # grids_list[[3]] = c(1:91)
    # grids_whole = expand.grid(grids_list)
    # grids = grids_whole[index,]

    n = 904
    V = nrow(grids)
    gibbs_c = chain$gibbs_c
    gibbs_tau_1_sq = chain$gibbs_tau_1_sq
    gibbs_tau_2_sq = chain$gibbs_tau_2_sq
    gibbs_thres = chain$gibbs_thres
    #gibbs_thres = rep(thres, T)
    temp_logL = chain$temp_logL
    rho_hat = matrix(nrow = V, ncol = T)
    rho_mean = rep(0, V)
    xi_hat_m = matrix(nrow = V, ncol = T)
    sigma_pos_sq_hat_m = matrix(nrow = V, ncol = T)
    sigma_neg_sq_hat_m = matrix(nrow = V, ncol = T)
    for(i in c(1:T)){
      xi_hat_m[,i] = Xmat %*% gibbs_c[,i]
      sigma_pos_sq_hat_m[,i] = ifelse(xi_hat_m[,i] > gibbs_thres[i],xi_hat_m[,i], 0)
      sigma_neg_sq_hat_m[,i] = ifelse(xi_hat_m[,i] < -gibbs_thres[i],-xi_hat_m[,i], 0)
    }
    #tau_1_sq = rowMeans(gibbs_tau_1_sq[,burn_in:T])
    #tau_2_sq = rowMeans(gibbs_tau_2_sq[,burn_in:T])
    log_tau_1_sq = Xmat%*%rnorm(length(lambda),mean=0,sd=0.5*sqrt(lambda))
    log_tau_2_sq = Xmat%*%rnorm(length(lambda),mean=0,sd=0.5*sqrt(lambda))
    tau_1_sq = exp(log_tau_1_sq)/5
    tau_2_sq = exp(log_tau_2_sq)/5
    c_init = rowMeans(chain$gibbs_c[,burn_in:T])
    sigma_pos_sq = rowMeans(sigma_pos_sq_hat_m[,burn_in:T])
    sigma_neg_sq = rowMeans(sigma_neg_sq_hat_m[,burn_in:T])

    sigma_pos = sqrt(sigma_pos_sq)
    sigma_neg = sqrt(sigma_neg_sq)
    eta_pos = matrix(NA, nrow=nrow(grids), ncol = n)
    eta_neg = matrix(NA, nrow=nrow(grids), ncol = n)
    Y_1 = matrix(NA, nrow=nrow(grids), ncol = n)
    Y_2 = matrix(NA, nrow=nrow(grids), ncol = n)
    for(i in 1:n){
      eta_pos[,i] = as.numeric(sigma_pos) * (Xmat %*% rnorm(length(lambda),mean = 0,sd = sqrt(lambda)))
      eta_neg[,i] = as.numeric(sigma_neg) * (Xmat %*% rnorm(length(lambda),mean = 0,sd = sqrt(lambda)))
      
    }
 
    for(i in 1:n){
      Y_1[,i] = rnorm(length(tau_1_sq), sd = sqrt(tau_1_sq))
      Y_2[,i] = rnorm(length(tau_2_sq), sd = sqrt(tau_2_sq))
    }
    Y_1 = Y_1 + eta_pos + eta_neg
    Y_2 = Y_2 + eta_pos - eta_neg
    thres_xi = sigma_pos_sq - sigma_neg_sq
    abs_thres_xi = abs(thres_xi)
    var_Y_1 = tau_1_sq + abs_thres_xi
    var_Y_2 = tau_2_sq + abs_thres_xi
    rho = thres_xi/sqrt(var_Y_1 * var_Y_2)
  
    cor_type = ifelse(rho > 0, 1, 0)
    cor_type = ifelse(rho < 0, -1, cor_type)
  

    ratio = (sigma_pos^2 + sigma_neg^2)/(sigma_pos^2 + sigma_neg^2 + tau_1_sq^2 + tau_2_sq^2) #### signal to noise ratio
    return(list('x' = grids, 'Xmat' = Xmat, 'lambda' = lambda, 'sigma_pos_sq' = sigma_pos_sq, 'sigma_neg_sq' = sigma_neg_sq, 'tau_1_sq' = tau_1_sq,
    'tau_2_sq' = tau_2_sq,'c_init' = c_init, 'rho' = rho, 'eta_pos' = eta_pos, 'eta_neg' = eta_neg,
     'Y_1' = Y_1, 'Y_2' = Y_2, 'var_Y_1' = var_Y_1, 'var_Y_2' = var_Y_2, 'ratio' = ratio, 'cor_type' = cor_type))
}
  

