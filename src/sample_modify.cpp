# include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
# include <math.h>
# include <list>
# include <iostream>
# include <vector>

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
mat fun_mul(mat A, colvec x) 
    { 
        A.each_col() %= x;
        return A;
    }
// [[Rcpp::export]]
mat fun_dev(mat A, colvec x) 
    { 
        A.each_col() /= x;
        return A;
    }
// [[Rcpp::export]]
vec G_thres_cpp(vec x, double thres){
	vec out(x.n_elem);
	for(int i=0; i < x.n_elem; i++){
		if(x(i) > thres){
			out(i) = x(i);
		}
		else{
			out(i) = 0;
		}
	}
	return(out);
}


// [[Rcpp::export]]

mat post_e_cpp(int l, mat Se_pos, mat Se_neg, mat Se_pos_other, mat Se_neg_other, vec Sc, vec r, vec s_sq, mat Xmat, vec lambda, mat Y_pos, mat Y_neg, double thres){
	mat res;
	mat res_pos;
	mat res_neg;
	vec Cl_pos = G_thres_cpp(Sc, thres) % Xmat.col(l);
	vec Cl_neg = G_thres_cpp(-Sc, thres) % Xmat.col(l);
	if(sum(Cl_pos == 0) == Cl_pos.size()){
		vec M_pos_0(Y_pos.n_cols, fill::zeros);
		vec Var_pos_0(Y_pos.n_cols, fill::zeros);
		res_pos = join_horiz(M_pos_0,Var_pos_0);
	}
	else{

		mat b_pos = fun_mul((Se_pos_other), G_thres_cpp(Sc, thres));
		mat mu_neg = fun_mul(Se_neg, G_thres_cpp(-Sc, thres));
		mat m_pos = fun_dev((Y_pos - b_pos) - fun_mul((Y_neg - mu_neg), r), Cl_pos);
		vec var_pos = (1.0 - square(r)) % s_sq / square(Cl_pos);
    	m_pos = m_pos.rows(find( Cl_pos != 0 ));
		var_pos = var_pos.elem(find(Cl_pos != 0));
		vec M_pos = sum((lambda(l) * fun_dev(m_pos, var_pos)).t(), 1)/as_scalar(1 + sum(lambda(l) / var_pos,0));
		vec Var_pos(M_pos.size(), 1, fill::ones);
		Var_pos = Var_pos * (lambda(l)/(1 + sum(lambda(l)/var_pos,0)));
		res_pos = join_horiz(M_pos,Var_pos);
	}
	if(sum(Cl_neg == 0) == Cl_neg.size()){
		vec M_neg_0(Y_neg.n_cols, fill::zeros);
		vec Var_neg_0(Y_neg.n_cols, fill::zeros);
		res_neg = join_horiz(M_neg_0,Var_neg_0);
	}
	else{
		mat b_neg = fun_mul((Se_neg_other), G_thres_cpp(-Sc, thres));
		mat mu_pos = fun_mul(Se_pos, G_thres_cpp(Sc, thres));
		mat m_neg = fun_dev((Y_neg - b_neg) - fun_mul((Y_pos - mu_pos), r), Cl_neg);
		vec var_neg = (1.0 - square(r)) % s_sq / square(Cl_neg);
    	m_neg = m_neg.rows(find( Cl_neg != 0 ));
		var_neg = var_neg.elem(find(Cl_neg != 0));
		vec M_neg = sum((lambda(l) * fun_dev(m_neg, var_neg)).t(), 1)/as_scalar(1 + sum(lambda(l) / var_neg,0));
		vec Var_neg(M_neg.size(), 1, fill::ones);
		Var_neg = Var_neg * (lambda(l)/(1 + sum(lambda(l)/var_neg,0)));
		res_neg = join_horiz(M_neg,Var_neg);
	}
	res =  join_horiz(res_pos, res_neg);
	return (res);
}

// [[Rcpp::export]]

mat post_tau_cpp(int n, mat Se_pos, mat Se_neg, vec Sc, mat Xmat, mat Y_pos, mat Y_neg, double alpha1, double beta1, double alpha2, double beta2, double thres){

	mat mu_pos = fun_mul(Se_pos, G_thres_cpp(Sc, thres));
	mat mu_neg = fun_mul(Se_neg, G_thres_cpp(-Sc, thres));
	mat T_pos = Y_pos - mu_pos;
	mat T_neg = Y_neg - mu_neg;
	vec shape1(Xmat.n_rows, 1, fill::ones);
	shape1 = shape1 * (alpha1 + double(n)/2);
	vec rate1 = sum(square(T_pos + T_neg),1)/2 + beta1;
	vec shape2(Xmat.n_rows, 1, fill::ones);
	shape2 = shape2 * (alpha2 + double(n)/2);
	vec rate2 = sum(square(T_pos - T_neg),1)/2 + beta2;
	mat res = join_horiz(shape1, rate1, shape2, rate2);
	return(res);
}


// [[Rcpp::export]]
double sample_c_l_cpp(mat grids, int n, int l, mat Se_pos, mat Se_neg, vec Sc_other, vec r, vec s_sq, double tau_xi_sq, mat Xmat, 
	vec lambda, mat Y_pos, mat Y_neg, double thres){
	//Rcout << "hello" << endl;

	int V = grids.n_rows;
	vec Lv = (thres - Sc_other) /  Xmat.col(l);
	vec Uv = (-thres - Sc_other) /  Xmat.col(l);
	vec Kv = (-2 * (1 - square(r))) % s_sq;
	vec Se_pos_sq = sum(square(Se_pos), 1);
	vec Se_neg_sq = sum(square(Se_neg), 1);
	vec Sc_sq = square(Sc_other);
	vec sign = Xmat.col(l);
	vec sign_2 = 2 * sign;
	vec sign_sq = square(sign);

	vec pos_pos = sum(Y_pos % Se_pos, 1);
	vec pos_neg = sum(Y_pos % Se_neg, 1);
	vec neg_pos = sum(Y_neg % Se_pos, 1);
	vec neg_neg = sum(Y_neg % Se_neg, 1);

	vec A_Pv = sum(Se_pos_sq, 1) % sign_sq;
	vec B_Pv = sign_2 % Sc_other % Se_pos_sq - sign_2 % pos_pos + (2 * r) % sign % neg_pos;
	vec C_Pv = Sc_sq % Se_pos_sq - (2 * pos_pos) % Sc_other + (2 * r) % neg_pos % Sc_other;

	vec A_Qv = sum(Se_neg_sq, 1) % sign_sq;
	vec B_Qv = sign_2 % Sc_other % Se_neg_sq + sign_2 % neg_neg - (2 * r) % sign % pos_neg;
	vec C_Qv = Sc_sq % Se_neg_sq + (2 * neg_neg) % Sc_other - (2 * r) % pos_neg % Sc_other;

	vec Lv_final = Lv.elem(find(sign !=0));
	vec Uv_final = Uv.elem(find(sign !=0));
	vec interval = sort(join_vert(Lv_final, Uv_final));
	interval.resize(interval.size() + 2);
	interval(interval.size() - 2) = min(interval) - 10000000;
	interval(interval.size() - 1) = max(interval) + 10000000;
	interval = sort(interval);
 	vec D(interval.size() - 1);
	vec E(interval.size() - 1);
	vec F(interval.size() - 1);
	vec M(interval.size() - 1);
	vec H(interval.size() - 1);
	double Dk; double Ek; double Fk;
	for(int k=0; k < (interval.size() - 1); k++){
		Dk = 0; Ek = 0; Fk = 0;
		for(int i = 0; i < V; i++){
			if(sign(i) > 0){ 
				if(interval(k) >= Lv(i)){
					Dk += (A_Pv(i)/double(Kv(i)));
					Ek += (B_Pv(i)/double(Kv(i)));
					Fk += (C_Pv(i)/double(Kv(i)));
				}
				if(interval(k+1) <= Uv(i)){
					Dk += (A_Qv(i)/double(Kv(i)));
					Ek += (B_Qv(i)/double(Kv(i)));
					Fk += (C_Qv(i)/double(Kv(i)));
				}
			}
			else if(sign(i) < 0){
				if(interval(k+1) <= Lv(i)){
					Dk += (A_Pv(i) /double(Kv(i)));
					Ek += (B_Pv(i) /double(Kv(i)));
					Fk += (C_Pv(i) /double(Kv(i)));
				}
				if(interval(k) >= Uv(i)){
					Dk += (A_Qv(i) /double(Kv(i)));
					Ek += (B_Qv(i) /double(Kv(i)));
					Fk += (C_Qv(i) /double(Kv(i)));
				}
			}
			else{
				
				continue;
			}
		}
		D(k) = Dk; E(k) = Ek; F(k) = Fk;
	}
	D = D - 1.0/(2 * tau_xi_sq * lambda(l));
	
	vec mu = - E / (2.0 * D);
	vec sig_sq = - 1 / (2.0 * D);

	for(int k = 0; k < (interval.size() - 1); k++){

		double La = R::pnorm((interval(k) - mu(k))/sqrt(sig_sq(k)), 0, 1, true, true);
		double Lb = R::pnorm((interval(k + 1) - mu(k))/sqrt(sig_sq(k)), 0, 1, true, true);
		H(k) = Lb + log(1-exp(La-Lb)) + F(k) - pow(E(k),2)/(4 * D(k)) + (double(1)/2) * log(double(1)/(-2 * tau_xi_sq * lambda(l) * D(k)));
		
	}

	double L0 = max(H);
	for(int k = 0; k < (interval.size() - 1); k++){
		M(k) = exp(H(k) - L0);
	}
	vec W = M/sum(M);
	double S = 0;
 	double rand = conv_to<double>::from(randu(1));
 	//Rcout << rand << endl;
 	//double rand = 0.5;
 	double res = 0;
 	for(int k=0; k < (interval.size() - 1); k++){
 		S += W(k);
 		if(W(k) == 0){
			continue;
		}
		if(k == 0){
			if(rand < S){
				res = as<double>(rtruncnorm(1, mu(k), sqrt(sig_sq(k)), interval(k), interval(k + 1)));
				//res = as<double>(rtruncnorm(1, interval(k), interval(k + 1), mu(k), sqrt(sig_sq(k))));
				break;
			}
		}
		else{
			if(rand > (S - W(k)) && rand < S){
				res = as<double>(rtruncnorm(1, mu(k), sqrt(sig_sq(k)), interval(k), interval(k + 1)));
				//res = as<double>(rtruncnorm(1, interval(k), interval(k + 1), mu(k), sqrt(sig_sq(k))));

				break;
			}
		}
 	}

 	return(res);
}




// [[Rcpp::export]]
double sample_thres_cpp(mat grids, int n, mat Se_pos, mat Se_neg, vec Sc, vec r, vec s_sq, 
	double tau_xi_sq, mat Xmat, vec lambda, mat Y_pos, mat Y_neg){

	int V = grids.n_rows;
	vec Sc_sq = square(Sc);
	vec Kv = (-2 * (1 - square(r))) % s_sq;
	//vec Pv = (2 * Sc) % (sum(Y_neg % Se_pos, 1) % r - sum(Y_pos % Se_pos, 1)) + square(Sc) % sum(square(Se_pos), 1);
	//vec Qv = (2 * Sc) % (sum(Y_neg % Se_neg, 1) - sum(Y_pos % Se_neg, 1) % r) + square(Sc) % sum(square(Se_neg), 1);
	vec Pv = (2 * Sc) % (fun_mul(sum(Y_neg % Se_pos, 1), r) - sum(Y_pos % Se_pos, 1)) + square(Sc) % sum(square(Se_pos), 1);
	vec Qv = (2 * Sc) % (sum(Y_neg % Se_neg, 1) - fun_mul(sum(Y_pos % Se_neg, 1), r)) + square(Sc) % sum(square(Se_neg), 1);
	vec interval = abs(Sc);
	interval = interval.elem(find(interval <= 3));
	interval.resize(interval.size() + 2);
	interval(interval.size() - 2) = 0;
	interval(interval.size() - 1) = 3;
	interval = sort(interval);
	vec D(interval.size() - 1);
	double Dk;
	for(int k = 0; k < interval.size() - 1; k++){
		Dk = 0;
		for(int i = 0; i < V; i++){
			if(interval(k + 1) <= Sc(i)){
				Dk = Dk + Pv(i)/Kv(i);
			}
			if(interval(k + 1) <= -Sc(i)){
 				Dk = Dk + Qv(i)/Kv(i);
 			}
		}
		D(k) = Dk;
 		D(k) = D(k) + log(interval(k+1) - interval(k)) + log(3);

	}
	double D0 = max(D);
	vec W = D - D0;
	W = exp(W);
	W = W/sum(W);
	double S = 0;
 	double rand = conv_to<double>::from(randu(1));
 	double res = 0;
 	for(int k=0; k < (interval.size() - 1); k++){
 		S += W(k);
 		if(W(k) == 0){
			continue;
		}
		if(k == 0){
			if(rand < S){
				res = (interval(k + 1) - interval(k)) * conv_to<double>::from(randu(1)) + interval(k);
				break;
			}
		}
		else{
			if(rand > (S - W(k)) && rand < S){
				res = (interval(k + 1) - interval(k)) * conv_to<double>::from(randu(1)) + interval(k);
				break;
			}
		}
 	}
 	return(res);	
}



// [[Rcpp::export]]

double loglikelihood_cpp(int n, mat Se_pos, mat Se_neg, vec Sc, vec r, vec s_sq, 
	mat Xmat, mat Y_pos, mat Y_neg, double thres){

	mat mu_pos = fun_mul(Se_pos, G_thres_cpp(Sc, thres));
	mat mu_neg = fun_mul(Se_neg, G_thres_cpp(-Sc, thres));
	vec H1 = sum(square(Y_pos - mu_pos) + square(Y_neg - mu_neg) - fun_mul((Y_pos - mu_pos) %  (Y_pos - mu_pos), (2 * r)), 1);
	vec H2 = 2 * (1 - square(r)) % s_sq;
	vec H3 = (2 * M_PI * s_sq) % sqrt(1 - square(r));
	vec logL = -n * log(H3) - H1/H2;
	return(sum(logL));
}


// [[Rcpp::export]]
List sample_gibbs_cpp(mat grids, int T, int V, int n, int L, mat Xmat, vec lambda,vec tau_1_sq_init, 
	vec tau_2_sq_init, vec c_init, mat e_pos_init, mat e_neg_init, mat Y_pos, mat Y_neg, double alpha1, 
	double beta1, double alpha2, double beta2, double thres_init, Function rinvgamma){
	//Rcout << "hello" << endl;
	mat temp_e_pos = e_pos_init;
	mat temp_e_neg = e_neg_init;
	vec temp_tau_1_sq = tau_1_sq_init;
	vec temp_tau_2_sq = tau_2_sq_init;
	vec temp_c = c_init;
	double temp_thres = thres_init;
	double temp_tau_xi_sq = 1;
	mat gibbs_tau_1_sq(V,T);
	mat gibbs_tau_2_sq(V,T);
	mat gibbs_c(L,T);
	gibbs_tau_1_sq.col(0) = temp_tau_1_sq;
	gibbs_tau_2_sq.col(0) = temp_tau_2_sq;

	gibbs_c.col(0) = temp_c;
	vector<double> temp_logL;
	vector<double> gibbs_thres;
	vec r = (temp_tau_1_sq - temp_tau_2_sq)/(temp_tau_1_sq + temp_tau_2_sq);
	vec s_sq = (temp_tau_1_sq + temp_tau_2_sq)/4;
	vec Sc = Xmat * temp_c;
	mat Se_pos = Xmat * temp_e_pos;
	mat Se_neg = Xmat * temp_e_neg;
	
	temp_logL.push_back(loglikelihood_cpp(n, Se_pos, Se_neg, Sc, r, s_sq, Xmat, Y_pos, Y_neg, temp_thres));
	gibbs_thres.push_back(thres_init);
	
	for(int t = 0; t < (T-1); t++){
		if(t % 1 == 0){
			Rcout << t << endl;
		}
		mat para_tau(V, 4);
		para_tau = post_tau_cpp(n, Se_pos, Se_neg, Sc, Xmat, Y_pos, Y_neg, alpha1, beta1, alpha2, beta2, temp_thres);
		/*gibbs_tau_1_sq.col(t+1) = randg(V) % (sqrt(para_tau.col(0)) % (1.0/para_tau.col(1))) + (1.0/para_tau.col(1)) % (para_tau.col(0) - sqrt(para_tau.col(0)));
		temp_tau_1_sq = gibbs_tau_1_sq.col(t+1);
		gibbs_tau_2_sq.col(t+1) = randg(V) % (sqrt(para_tau.col(2)) % (1.0/para_tau.col(3))) + (1.0/para_tau.col(3)) % (para_tau.col(2) - sqrt(para_tau.col(2)));
		temp_tau_2_sq = gibbs_tau_2_sq.col(t+1);*/
		/*for(int v=0; v < V; v++){

			gibbs_tau_1_sq(v, t+1) = 1.0/randg<double>(distr_param(para_tau(v, 0), 1.0/para_tau(v, 1)));
			temp_tau_1_sq(v) = gibbs_tau_1_sq(v, t+1);
			gibbs_tau_2_sq(v, t+1) = 1.0/randg<double>(distr_param(para_tau(v, 2), 1.0/para_tau(v, 3)));
			temp_tau_2_sq(v) = gibbs_tau_2_sq(v, t+1);
		}*/
		gibbs_tau_1_sq.col(t+1) = as<vec>(rinvgamma(V, para_tau.col(0), para_tau.col(1)));
		temp_tau_1_sq = gibbs_tau_1_sq.col(t+1);
		gibbs_tau_2_sq.col(t+1) = as<vec>(rinvgamma(V, para_tau.col(2), para_tau.col(3)));
		temp_tau_2_sq = gibbs_tau_2_sq.col(t+1);
		r = (temp_tau_1_sq - temp_tau_2_sq)/(temp_tau_1_sq + temp_tau_2_sq);
		s_sq = (temp_tau_1_sq + temp_tau_2_sq)/4;
		for(int l = 0; l < L; l++){
			// sample c
			vec Sc_other = Sc - temp_c(l) * Xmat.col(l);
			gibbs_c(l,t+1) = sample_c_l_cpp(grids, n, l, Se_pos, Se_neg, Sc_other, r, s_sq, temp_tau_xi_sq, Xmat, lambda, Y_pos, Y_neg, temp_thres);
			temp_c(l) = gibbs_c(l,t+1);
			Sc = Sc - (gibbs_c(l,t) - gibbs_c(l,t+1)) * Xmat.col(l);
			// sample e
			mat para_e(n, 4);
			mat Se_pos_other = Se_pos - Xmat.col(l) * temp_e_pos.row(l);
			mat Se_neg_other = Se_neg - Xmat.col(l) * temp_e_neg.row(l);
			para_e = post_e_cpp(l, Se_pos, Se_neg,  Se_pos_other, Se_neg_other, Sc, r, s_sq, Xmat, lambda, Y_pos, Y_neg, temp_thres);
			mat temp_e_pos_before_update = temp_e_pos.row(l);
			mat temp_e_neg_before_update = temp_e_neg.row(l);
			temp_e_pos.row(l) = (randn(n) % sqrt(para_e.col(1)) + para_e.col(0)).as_row();
			temp_e_neg.row(l) = (randn(n) % sqrt(para_e.col(3)) + para_e.col(2)).as_row();
			//Rcout << l << endl;
			Se_pos = Se_pos + Xmat.col(l) * (temp_e_pos.row(l) - temp_e_pos_before_update);
			Se_neg = Se_neg + Xmat.col(l) * (temp_e_neg.row(l) - temp_e_neg_before_update);
		}
		
		temp_thres = sample_thres_cpp(grids, n, Se_pos, Se_neg, Sc, r, s_sq, temp_tau_xi_sq, Xmat, lambda, Y_pos, Y_neg);
		gibbs_thres.push_back(temp_thres);
		temp_logL.push_back(loglikelihood_cpp(n, Se_pos, Se_neg, Sc, r, s_sq, Xmat, Y_pos, Y_neg, temp_thres));
	}
	List chain = List::create(Named("gibbs_c") = gibbs_c , _["gibbs_tau_1_sq"] = gibbs_tau_1_sq, _["gibbs_tau_2_sq"] = gibbs_tau_2_sq, _["temp_logL"] = temp_logL, _["gibbs_thres"] = gibbs_thres);
	//List chain = List::create(Named("gibbs_c") = gibbs_c , _["gibbs_tau_1_sq"] = gibbs_tau_1_sq, _["gibbs_tau_2_sq"] = gibbs_tau_2_sq);

	return (chain);
}




/*
post_tau_cpp(100, e_pos = dat$true_e_pos, e_neg = dat$true_e_neg, dat$true_c, Xmat, Y_pos, Y_neg, 1, 1, 1, 1, 1)
*/
