// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace arma;
using namespace std;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
vec rowsum_Mat(mat M) {
    int nr=M.n_rows;
    vec out(nr);
    for(int i=0;i<nr;i++){
        out(i)=sum(M.row(i));
    }
    return out;
}

// [[Rcpp::export]]
vec colsum_Mat(mat M) {
    int nc=M.n_cols;
    vec out(nc);
    for(int i=0;i<nc;i++){
        out(i)=sum(M.col(i));
    }
    return out;
}

// [[Rcpp::export]]
float epan(float input){
    float output;
    if(abs(input)<=1){
        output=0.75*(1-pow(input,2));
    }
    else{
        output=0;
    }
    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// gamma, Tau update function, gradient and hessian functions and ELBO convergence function

// [[Rcpp::export]]
cube gamma_update_HMM_stat_dir(mat gamma, vec pi, mat theta, mat network, int N, int K){
    cube quad_lin_coeff(N,K,2);
    for(int i = 0; i < N; i++){
        if(i!=(N-1)){
            for(int k = 0; k < K; k++){
                float t1=0;
                for(int j = i+1; j < N; j++){
                    for(int l = 0; l < K; l++){
                        float exp_val_1=exp(theta(k,0));
                        float exp_val_2=exp(theta(l,0));
                        float exp_val_3=exp(theta(k,1)+theta(l,1));
                        vec alpha(4,fill::zeros);
                        alpha(1)=exp_val_1;
                        alpha(2)=exp_val_2;
                        alpha(3)=exp_val_3;
                        float alpha_max=alpha.max();
                        float exp_val_1_mod=exp(theta(k,0)-alpha_max);
                        float exp_val_2_mod=exp(theta(l,0)-alpha_max);
                        float exp_val_3_mod=exp(theta(k,1)+theta(l,1)-alpha_max);
                        float log_exp_val=alpha_max+log(exp(-alpha_max)+exp_val_1_mod+exp_val_2_mod+exp_val_3_mod);
                        int indicator_10=(network(i,j)==1)&(network(j,i)==0);
                        int indicator_01=(network(i,j)==0)&(network(j,i)==1);
                        int indicator_11=(network(i,j)==1)&(network(j,i)==1);
                        t1+=((gamma(j,l)/(2*gamma(i,k)))*((indicator_10*theta(k,0))+(indicator_01*theta(l,0))+((indicator_11)*(theta(k,1)+theta(l,1)))-log_exp_val));
                    }
                }
                quad_lin_coeff(i,k,0)=t1-(1/gamma(i,k));
                quad_lin_coeff(i,k,1)=(log(pi(k))-log(gamma(i,k))+1);
            }
        } else if(i==(N-1)){
            for(int k = 0; k < K; k++){
                quad_lin_coeff(i,k,0)=-(1/gamma((N-1),k));
                quad_lin_coeff(i,k,1)=(log(pi(k))-log(gamma((N-1),k))+1);
            }
        }
    }
    return quad_lin_coeff;
}

// [[Rcpp::export]]
vec grad_HMM_stat_dir_oe(mat theta, mat gamma, mat network, int N, int K){
    vec grad_vector(K,fill::zeros);
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            mat grad_mat_1(K,K,fill::zeros);
            mat grad_mat_2(K,K,fill::zeros);
            for(int k = 0; k < K; k++){
                for(int l = 0; l < K; l++){
                    float exp_val_1=exp(theta(k,0));
                    float exp_val_2=exp(theta(l,0));
                    float exp_val_3=exp(theta(k,1)+theta(l,1));
                    int indicator_10=(network(i,j)==1)&(network(j,i)==0);
                    int indicator_01=(network(i,j)==0)&(network(j,i)==1);
                    grad_mat_1(k,l)=gamma(i,k)*gamma(j,l)*(indicator_10-(exp_val_1/(1+exp_val_1+exp_val_2+exp_val_3)));
                    grad_mat_2(k,l)=gamma(i,k)*gamma(j,l)*(indicator_01-(exp_val_2/(1+exp_val_1+exp_val_2+exp_val_3)));
                }
            }
            vec rsum_mat_1=rowsum_Mat(grad_mat_1);
            vec csum_mat_2=colsum_Mat(grad_mat_2);
            grad_vector+=rsum_mat_1+csum_mat_2;
        }
    }
    return grad_vector;
}

// [[Rcpp::export]]
vec grad_HMM_stat_dir_re(mat theta, mat gamma, mat network, int N, int K){
    vec grad_vector(K,fill::zeros);
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            mat grad_mat(K,K,fill::zeros);
            for(int k = 0; k < K; k++){
                for(int l = 0; l < K; l++){
                    float exp_val_1=exp(theta(k,0));
                    float exp_val_2=exp(theta(l,0));
                    float exp_val_3=exp(theta(k,1)+theta(l,1));
                    int indicator_11=(network(i,j)==1)&(network(j,i)==1);
                    grad_mat(k,l)=gamma(i,k)*gamma(j,l)*(indicator_11-(exp_val_3/(1+exp_val_1+exp_val_2+exp_val_3)));
                }
            }
            vec rsum_mat=rowsum_Mat(grad_mat);
            vec csum_mat=colsum_Mat(grad_mat);
            grad_vector+=rsum_mat+csum_mat;
        }
    }
    return grad_vector;
}

// [[Rcpp::export]]
mat hess_HMM_stat_dir_oe(mat theta, mat gamma, int N, int K){
    mat t1(K,K);
    mat hess_mat_1(K,K,fill::zeros);
    mat hess_mat_2(K,K,fill::zeros);
    mat hess_mat_3(K,K,fill::zeros);
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            mat hess_matsub_1(K,K);
            mat hess_matsub_2(K,K);
            mat hess_matsub_3(K,K);
            for(int k = 0; k < K; k++){
                for(int l = 0; l < K; l++){
                    float exp_val_1=exp(theta(k,0));
                    float exp_val_2=exp(theta(l,0));
                    float exp_val_3=exp(theta(k,1)+theta(l,1));
                    hess_matsub_1(k,l)=(gamma(i,k)*gamma(j,l)*((exp_val_1*exp_val_2)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2))));
                    hess_matsub_2(k,l)=(gamma(i,k)*gamma(j,l)*((exp_val_1+exp_val_1*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2))));
                    hess_matsub_3(k,l)=(gamma(i,k)*gamma(j,l)*((exp_val_2+exp_val_2*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2))));
                }
            }
            hess_mat_1+=hess_matsub_1;
            hess_mat_2+=hess_matsub_2;
            hess_mat_3+=hess_matsub_3;
        }
    }
    for(int k = 0; k < K; k++){
        for(int l = 0; l < K; l++){
            if(k!=l){
                t1(k,l)=(hess_mat_1(k,l)+hess_mat_1(l,k));
            }
        }
    }
    vec rsum_1=rowsum_Mat(hess_mat_1);
    vec csum_1=colsum_Mat(hess_mat_1);
    vec rsum_2=rowsum_Mat(hess_mat_2);
    vec csum_3=colsum_Mat(hess_mat_3);
    for(int k = 0; k < K; k++){
        t1(k,k)=(-csum_3(k)-rsum_2(k)-(rsum_1(k)+csum_1(k)-2*hess_mat_1(k,k)));
    }
    return t1;
}

// [[Rcpp::export]]
mat hess_HMM_stat_dir_re(mat theta, mat gamma, int N, int K){
    mat t1(K,K);
    mat hess_mat(K,K,fill::zeros);
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            mat hess_matsub(K,K);
            for(int k = 0; k < K; k++){
                for(int l = 0; l < K; l++){
                    float exp_val_1=exp(theta(k,0));
                    float exp_val_2=exp(theta(l,0));
                    float exp_val_3=exp(theta(k,1)+theta(l,1));
                    hess_matsub(k,l)=(gamma(i,k)*gamma(j,l)*(((1+exp_val_1+exp_val_2)*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2))));
                }
            }
            hess_mat+=hess_matsub;
        }
    }
    for(int k = 0; k < K; k++){
        for(int l = 0; l < K; l++){
            if(k!=l){
                t1(k,l)=-(hess_mat(k,l)+hess_mat(l,k));
            }
        }
    }
    vec rsum=rowsum_Mat(hess_mat);
    vec csum=colsum_Mat(hess_mat);
    for(int k = 0; k < K; k++){
        t1(k,k)=(-csum(k)-rsum(k)-2*hess_mat(k,k));
    }
    return t1;
}

// [[Rcpp::export]]
mat hess_HMM_stat_dir_oe_re(mat theta, mat gamma, int N, int K){
    mat t1(K,K);
    mat hess_mat_1(K,K,fill::zeros);
    mat hess_mat_2(K,K,fill::zeros);
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            mat hess_matsub_1(K,K);
            mat hess_matsub_2(K,K);
            for(int k = 0; k < K; k++){
                for(int l = 0; l < K; l++){
                    float exp_val_1=exp(theta(k,0));
                    float exp_val_2=exp(theta(l,0));
                    float exp_val_3=exp(theta(k,1)+theta(l,1));
                    hess_matsub_1(k,l)=(gamma(i,k)*gamma(j,l)*((exp_val_1*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2))));
                    hess_matsub_2(k,l)=(gamma(i,k)*gamma(j,l)*((exp_val_2*exp_val_3)/(pow((1+exp_val_1+exp_val_2+exp_val_3),2))));
                }
            }
            hess_mat_1+=hess_matsub_1;
            hess_mat_2+=hess_matsub_2;
        }
    }
    for(int k = 0; k < K; k++){
        for(int l = 0; l < K; l++){
            if(k!=l){
                t1(k,l)=(hess_mat_1(k,l)+hess_mat_2(l,k));
            }
        }
    }
    vec rsum_1=rowsum_Mat(hess_mat_1);
    vec csum_2=colsum_Mat(hess_mat_2);
    for(int k = 0; k < K; k++){
        t1(k,k)=(csum_2(k)+rsum_1(k)+(hess_mat_1(k,k)+hess_mat_2(k,k)));
    }
    return t1;
}

// [[Rcpp::export]]
float ELBO_conv_HMM_stat_dir(mat gamma, vec alpha, mat theta, mat network, int N, int K){
    float t1=0;
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            for(int k = 0; k < K; k++){
                for(int l = 0; l < K; l++){
                    float exp_val_1=exp(theta(k,0));
                    float exp_val_2=exp(theta(l,0));
                    float exp_val_3=exp(theta(k,1)+theta(l,1));
                    int indicator_10=(network(i,j)==1)&(network(j,i)==0);
                    int indicator_01=(network(i,j)==0)&(network(j,i)==1);
                    int indicator_11=(network(i,j)==1)&(network(j,i)==1);
                    t1+=(gamma(i,k)*gamma(j,l)*((indicator_10*theta(k,0))+(indicator_01*theta(l,0))+(indicator_11*(theta(k,1)+theta(l,1)))-log(1+exp_val_1+exp_val_2+exp_val_3)));
                }
            }
        }
    }
    
    
    float t2=0;
    for(int i = 0; i < N; i++){
        for(int k = 0; k < K; k++){
            if((alpha(k)>=(pow(10,(-100))))&(gamma(i,k)>=(pow(10,(-100))))){
                t2+=gamma(i,k)*(log(alpha(k))-log(gamma(i,k)));
            }
        }
    }
    
    float ELBO_val=t1+t2;
    return ELBO_val;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// Defining functions for K=1

// [[Rcpp::export]]
float grad_HMM_stat_dir_oe_K1(vec theta, mat network, int N){
    float exp_val_1=exp(theta(0));
    float exp_val_2=exp(2*theta(1));
    float exp_val_ratio=(2*exp_val_1)/(1+2*exp_val_1+exp_val_2);
    float grad_val_oe=0;
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            int indicator_10=(network(i,j)==1)&(network(j,i)==0);
            int indicator_01=(network(i,j)==0)&(network(j,i)==1);
            grad_val_oe+=((indicator_10+indicator_01)-exp_val_ratio);
        }
    }
    return grad_val_oe;
}

// [[Rcpp::export]]
float grad_HMM_stat_dir_re_K1(vec theta, mat network, int N){
    float exp_val_1=exp(theta(0));
    float exp_val_2=exp(2*theta(1));
    float exp_val_ratio=(2*exp_val_2)/(1+2*exp_val_1+exp_val_2);
    float grad_val_re=0;
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            int indicator_11=(network(i,j)==1)&(network(j,i)==1);
            grad_val_re+=((2*indicator_11)-exp_val_ratio);
        }
    }
    return grad_val_re;
}

// [[Rcpp::export]]
float hess_HMM_stat_dir_oe_K1(vec theta, int N){
    float exp_val_1=exp(theta(0));
    float exp_val_2=exp(2*theta(1));
    float exp_val_ratio=(exp_val_1+exp_val_1*exp_val_2)/(pow((1+2*exp_val_1+exp_val_2),2));
    float hess_val=0;
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            hess_val+=exp_val_ratio;
        }
    }
    float hess_val_oe=-2*hess_val;
    return hess_val_oe;
}

// [[Rcpp::export]]
float hess_HMM_stat_dir_re_K1(vec theta, int N){
    float exp_val_1=exp(theta(0));
    float exp_val_2=exp(2*theta(1));
    float Num=4*(exp_val_2+(2*exp_val_1*exp_val_2));
    float Denom=pow((1+2*exp_val_1+exp_val_2),2);
    float exp_val_ratio=Num/Denom;
    float hess_val_re=0;
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            hess_val_re+=(-exp_val_ratio);
        }
    }
    return hess_val_re;
}

// [[Rcpp::export]]
float hess_HMM_stat_dir_oe_re_K1(vec theta, int N){
    float exp_val_1=exp(theta(0));
    float exp_val_2=exp(2*theta(1));
    float exp_val_ratio=(4*exp_val_1*exp_val_2)/(pow((1+2*exp_val_1+exp_val_2),2));
    float hess_val_oe_re=0;
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            hess_val_oe_re+=exp_val_ratio;
        }
    }
    return hess_val_oe_re;
}

// [[Rcpp::export]]
float ELBO_conv_HMM_stat_dir_K1(vec theta, mat network, int N){
    float exp_val_1=exp(theta(0));
    float exp_val_2=exp(2*theta(1));
    float log_exp_val=log(1+2*exp_val_1+exp_val_2);
    float ELBO_val=0;
    for(int i = 0; i < (N-1); i++){
        for(int j = i+1; j < N; j++){
            int indicator_10=(network(i,j)==1)&(network(j,i)==0);
            int indicator_01=(network(i,j)==0)&(network(j,i)==1);
            int indicator_11=(network(i,j)==1)&(network(j,i)==1);
            ELBO_val+=((indicator_10*theta(0))+(indicator_01*theta(0))+(indicator_11*(2*theta(1)))-log_exp_val);
        }
    }
    return ELBO_val;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
