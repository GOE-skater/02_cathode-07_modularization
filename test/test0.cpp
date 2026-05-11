//*****************************************************************
//*****************************************************************
//**                                                             **
//**          plasma_check                                       **
//**          coded by Ryo Shirakawa.                            **
//*****************************************************************
//*****************************************************************

#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include "common.hpp"
#include "functions.hpp"
#include <sstream>

//*****************************************************************
//**                                                             **
//**           void solve_rhoe_simple                         **
//**                                                             **
//*****************************************************************
void solve_rhoe_wdTe_wSEE_PC(){

    std::vector<std::vector<double> > aW(ni+2,std::vector<double>(nj+2,0.0)); //係数
    std::vector<std::vector<double> > aE(ni+2,std::vector<double>(nj+2,0.0)); //係数
    std::vector<std::vector<double> > aS(ni+2,std::vector<double>(nj+2,0.0)); //係数
    std::vector<std::vector<double> > aN(ni+2,std::vector<double>(nj+2,0.0)); //係数
    std::vector<std::vector<double> > aP(ni+2,std::vector<double>(nj+2,0.0)); //係数
    std::vector<std::vector<double> > b (ni+2,std::vector<double>(nj+2,0.0)); //ソース
    std::vector<std::vector<double> > rhoe_new(ni+2,std::vector<double>(nj+2,0.0)); //ソース
    std::vector<std::vector<double> > rhoe_old_tmp(ni+2,std::vector<double>(nj+2,0.0));
    std::vector<std::vector<double> > rhoUex_old_tmp(ni+2,std::vector<double>(nj+2,0.0));
    std::vector<std::vector<double> > rhoUer_old_tmp(ni+2,std::vector<double>(nj+2,0.0));
    std::vector<std::vector<double> > drhoe(ni+2,std::vector<double>(nj+2,0.0));

    //Numerical treatment for exponential overflow and zero-division
    //------------------------------------
    double lim_exp = 1.0e2; //expの制限
    double EPS = 1.0e-3; //ゼロ割回避用の数字
    //------------------------------------

    //Acceleration factor for PC method
    //------------------------------------
    double c_mu = 1.0;
    //------------------------------------

    //convergence criteria
    //------------------------------------
    int icon_error2 = 1; //擬似ステップのエラー収束判定 0:max, 1:RMS
    //------------------------------------

    //output duration
    //------------------------------------
    int ndiv_out = 100;
    //------------------------------------
    
    //Reserve the past values
    //------------------------------------
    for (int i=0;i<=ni+1;i++){ 
        for (int j=0;j<=nj+1;j++){ 
            rhoe_old2[i][j] = rhoe_old[i][j]; //電場スキーム確認用
            rhoe_old[i][j] = rhoe[i][j];
            rhoUex_old[i][j] = rhoUex[i][j];
            rhoUer_old[i][j] = rhoUer[i][j];
        }
    }
    //------------------------------------

    //Initial prediction of correction
    //------------------------------------
    for (int i=0;i<=ni+1;i++){ 
        for (int j=0;j<=nj+1;j++){ 
            drhoe[i][j] = rhoe[i][j] - rhoe_old[i][j];
        }
    }
    //------------------------------------

    //calculate maximum electron density to avoid zero-division in Te calc.
    //------------------------------------
    double rhoe_max = 0.0;
    for (int iblock=0;iblock<2;iblock++){ 
        for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
            for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                if(rhoe_max < rhoe[i][j]){
                    rhoe_max = rhoe[i][j];
                }
            }
        }
    }
    //------------------------------------

    //PC method iteration
    //------------------------------------
    int ncount = 0; //number of iteration of PC method
    int icon_end = 0;
    do{
        
        ncount ++;

        //Reserve previous-step values
        //------------------------------------
        for (int i=0;i<=ni+1;i++){ 
            for (int j=0;j<=nj+1;j++){ 
                rhoe_old_tmp[i][j] = rhoe[i][j];
                rhoUex_old_tmp[i][j] = rhoUex[i][j];
                rhoUer_old_tmp[i][j] = rhoUer[i][j];
            }
        }
        //------------------------------------

        //Predict flux
        //==============================================
        //for rhoUex
        //------------------------------------
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flx_bl[iblock][0];i<=i_flx_bl[iblock][1];i++){
                for (int j=j_flx_bl[iblock][0];j<=j_flx_bl[iblock][1];j++){
                    
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rL = (r[j]+r[j-1])/2.0;
                    double rR = (r[j]+r[j+1])/2.0;
                    //------------------------------------

                    //moblity at cell interface
                    //------------------------------------
                    double mu_para_tmp = (mu_para[i][j] + mu_para[i-1][j])/2.0;
                    double mu_perp_tmp = (mu_perp[i][j] + mu_perp[i-1][j])/2.0;
                    
                    double mu_star_tmp = mu_perp_tmp*c_mu;
                    //------------------------------------

                    //calculate inverse moblity tensor at cell interface
                    //------------------------------------
                    double mu_rr_tmp = (mu_rr[i][j] + mu_rr[i-1][j])/2.0;
                    double mu_xr_tmp = (mu_xr[i][j] + mu_xr[i-1][j])/2.0;

                    double mu_inv_xx_tmp = mu_rr_tmp/(mu_para_tmp*mu_perp_tmp);
                    double mu_inv_xr_tmp = -mu_xr_tmp/(mu_para_tmp*mu_perp_tmp);

                    //double mu_xx_tmp = (mu_xx[i][j] + mu_xx[i-1][j])/2.0;
                    //double mu_rr_tmp = (mu_rr[i][j] + mu_rr[i-1][j])/2.0;
                    //double mu_xr_tmp = (mu_xr[i][j] + mu_xr[i-1][j])/2.0;

                    //double deno = mu_xx_tmp*mu_rr_tmp - mu_xr_tmp*mu_xr_tmp;
                    //double mu_inv_xx_tmp = mu_rr_tmp/(deno + 1e-100);
                    //double mu_inv_xr_tmp = -mu_xr_tmp/(deno + 1e-100);
                    //------------------------------------

                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double Te_tmp = (Te[i][j] + Te[i-1][j])/2.0; //electron temperature
                    double dTedx_tmp = (Te[i][j] - Te[i-1][j])/dx*double(icon_dTe); //electron temperature gradient
                    double vd_tmp = mu_star_tmp*(-Ex[i][j] - Boltz/e0*dTedx_tmp) + EPS; //electron drift velocity
                    double D_tmp = Boltz*Te_tmp/e0*mu_star_tmp; //electron diffusion coefficient
                    double z_tmp = vd_tmp*dx/D_tmp; //Peclet number
                    //------------------------------------

                    //interpolation of transverse flux component (should be old value)
                    //------------------------------------
                    double rhoUer_tmp = (rR*rhoUer_old_tmp[i][j+1] + rL*rhoUer_old_tmp[i][j]
                        + rR*rhoUer_old_tmp[i-1][j+1] + rL*rhoUer_old_tmp[i-1][j])/(4.0*r[j]);
                    //------------------------------------

                    //update flux
                    //------------------------------------
                    rhoUex[i][j] =
                          vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*rhoe[i-1][j]
                        - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*rhoe[i][j] 
                        + rhoUex_old_tmp[i][j] 
                        - mu_star_tmp*(mu_inv_xx_tmp*rhoUex_old_tmp[i][j] + mu_inv_xr_tmp*rhoUer_tmp);
                    //------------------------------------
                }
            }     
        }
        //------------------------------------

        //left-wall-flux (chamber)
        //------------------------------------
        for (int j=j_flc_bl[0][0];j<=j_flc_bl[0][1];j++){
            int i=i_flc_bl[0][0];

            //bulk-flux
            //------------------------------------
            rhoUex[i][j] = 0.0;
            //------------------------------------

            //wall-flux
            //------------------------------------
            double vth = sqrt(8.0*Boltz*Te[i][j]/(M_PI*masse));
            double delta_phi = fmax(-Ex[i][j]*dx/2.0,0.0); //adjacent cell-center potential from wall
            double coef_angle = 1.0;

            if(icon_eKineBC == 2) delta_phi = 0.0; //no need to extrapolate density
            if(icon_eFLRBC == 1){
                double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
                double angle = asin(fabs(Bx[i][j]/Bmag)); //B-field angle from x-axis
                double vth2 = sqrt(2.0*Boltz*Te[i][j]/masse);
                double omega_ce = e0*Bmag/masse;
                double local_r = vth2/omega_ce;
                coef_angle = sqrt(pow(sin(angle),2) + pow(fmin(local_r/dx,1.0)*cos(angle),2));
            }else if(icon_eFLRBC == 2){
                double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
                double angle = asin(fabs(Bx[i][j]/Bmag)); //B-field angle from x-axis
                double Halle = e0*Bmag/(masse*nu_m[i][j]);
                coef_angle = pow(sin(angle),2) + pow(cos(angle),2)/(1.0 + Halle*Halle);
            }

            double GammaPerN_wall_Lx = 
                -0.25*vth*coef_angle
                *exp(fmin(-e0*delta_phi/(Boltz*Te[i][j] + 1e-100),lim_exp))
                *(exp(-Utex[i][j]*Utex[i][j]) + utex[i][j]*sqrt(M_PI)*(erf(Utex[i][j]) - 1.0));

            rhoUex_wall[i][j] = rhoe[i][j]*GammaPerN_wall_Lx;
            //------------------------------------
        }
        //------------------------------------

        //left-wall-flux (magnet)
        //------------------------------------
        for (int j=j_flc_bl[1][0];j<=j_flc_bl[0][0]-1;j++){
            int i=i_flc_bl[1][0];

            //bulk-flux
            //------------------------------------
            rhoUex[i][j] = 0.0;
            //------------------------------------

            //wall-flux
            //------------------------------------
            double vth = sqrt(8.0*Boltz*Te[i][j]/(M_PI*masse));
            double delta_phi = fmax(-Ex[i][j]*dx/2.0,0.0); //adjacent cell-center potential from wall
            double coef_angle = 1.0;
            if(icon_eKineBC == 2) delta_phi = 0.0; //no need to extrapolate density
            if(icon_eFLRBC == 1){
                double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
                double angle = asin(fabs(Bx[i][j]/Bmag)); //B-field angle from x-axis
                double vth2 = sqrt(2.0*Boltz*Te[i][j]/masse);
                double omega_ce = e0*Bmag/masse;
                double local_r = vth2/omega_ce;
                coef_angle = sqrt(pow(sin(angle),2) + pow(fmin(local_r/dx,1.0)*cos(angle),2));
            }else if(icon_eFLRBC == 2){
                double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
                double angle = asin(fabs(Bx[i][j]/Bmag)); //B-field angle from x-axis
                double Halle = e0*Bmag/(masse*nu_m[i][j]);
                coef_angle = pow(sin(angle),2) + pow(cos(angle),2)/(1.0 + Halle*Halle);
            }
            double GammaPerN_wall_Lx = 
                -0.25*vth*coef_angle
                *exp(fmin(-e0*delta_phi/(Boltz*Te[i][j] + 1e-100),lim_exp))
                *(exp(-Utex[i][j]*Utex[i][j]) + utex[i][j]*sqrt(M_PI)*(erf(Utex[i][j]) - 1.0));

            rhoUex_wall[i][j] = rhoe[i][j]*GammaPerN_wall_Lx;
            //------------------------------------
        }
        //------------------------------------

        if(icon_open == 0){
            //right-wall-flux
            //------------------------------------
            for (int j=j_flc_bl[1][0];j<=j_flc_bl[1][1];j++){
                int i=i_flc_bl[1][1];

                //bulk-flux
                //------------------------------------
                rhoUex[i+1][j] = 0.0;
                //------------------------------------

                //calulate kinetic vector flux
                //------------------------------------
                double vth = sqrt(8.0*Boltz*Te[i][j]/(M_PI*masse));
                double delta_phi =  fmax(Ex[i+1][j]*dx/2.0,0.0); //adjacent cell-center potential from wall
                double coef_angle = 1.0;
                if(icon_eKineBC == 2) delta_phi = 0.0; //no need to extrapolate density
                if(icon_eFLRBC == 1){
                    double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
                    double angle = asin(fabs(Bx[i][j]/Bmag)); //B-field angle from x-axis
                    double vth2 = sqrt(2.0*Boltz*Te[i][j]/masse);
                    double omega_ce = e0*Bmag/masse;
                    double local_r = vth2/omega_ce;
                    coef_angle = sqrt(pow(sin(angle),2) + pow(fmin(local_r/dx,1.0)*cos(angle),2));
                }else if(icon_eFLRBC == 2){
                    double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
                    double angle = asin(fabs(Bx[i][j]/Bmag)); //B-field angle from x-axis
                    double Halle = e0*Bmag/(masse*nu_m[i][j]);
                    coef_angle = pow(sin(angle),2) + pow(cos(angle),2)/(1.0 + Halle*Halle);
                }
                double GammaPerN_wall_Rx = 
                    0.25*vth*coef_angle
                    *exp(fmin(-e0*delta_phi/(Boltz*Te[i][j] + 1e-100),lim_exp))
                    *(exp(-Utex[i+1][j]*Utex[i+1][j]) + utex[i+1][j]*sqrt(M_PI)*(erf(Utex[i+1][j]) + 1.0));

                rhoUex_wall[i+1][j] =  rhoe[i][j]*GammaPerN_wall_Rx;
                //------------------------------------
            }
            //------------------------------------

        }else{
            if(icon_opGauss == 0){
                //right-open-flux
                //------------------------------------
                for (int j=j_flc_bl[1][0];j<=j_flc_bl[1][1];j++){
                    int i=i_flc_bl[1][1];

                    //bulk-flux (equal to ion)
                    //------------------------------------
                    double rhoUix_Rx = rhoi[i][j]*Uix[i+1][j];
                    if(icon_ionConsv == 1){
                        rhoUix_Rx = rhoUix[i+1][j];
                    }
                    double Gamma_open_Rx = rhoUix_Rx;
                    rhoUex[i+1][j] =  Gamma_open_Rx;
                    //------------------------------------

                    //wall-flux
                    //------------------------------------
                    rhoUex_wall[i+1][j] = 0.0;
                    //------------------------------------
                }
                //------------------------------------
            }else{
                //right-open-flux
                //------------------------------------
                for (int j=j_flc_bl[1][0];j<=j_flc_bl[1][1];j++){
                    int i=i_flc_bl[1][1];

                    //outgoing flux
                    //------------------------------------
                    double vth = sqrt(8.0*Boltz*Te[i][j]/(M_PI*masse));
                    double GammaPerN_wall_Rx = 
                        0.25*vth
                        *(exp(-utex[i+1][j]*utex[i+1][j]) + utex[i+1][j]*sqrt(M_PI)*(erf(utex[i+1][j]) + 1.0));

                    double rhoUex_out = rhoe[i][j]*GammaPerN_wall_Rx;
                    //------------------------------------

                    //ion bulk-flux
                    //------------------------------------
                    double rhoUix_Rx = rhoi[i][j]*Uix[i+1][j];
                    if(icon_ionConsv == 1){
                        rhoUix_Rx = rhoUix[i+1][j];
                    }
                    //------------------------------------

                    //incoming flux
                    //------------------------------------
                    double rhoUex_in = fmin(rhoUix_Rx - eps0/(e0*dt)*(1.0 - R_Eop)*Ex_old[i+1][j] - rhoUex_out,0.0);
                    //std::cout << itime << ","<< rhoUix_Rx - eps0/(e0*dt)*(1.0 - R_Eop)*Ex_old[i+1][j] - rhoUex_out << std::endl;
                    //double rhoUex_in = rhoUix_Rx - eps0/(e0*dt)*(1.0 - R_Eop)*Ex_old[i+1][j] - rhoUex_out;
                    //------------------------------------

                    //bulk-flux
                    //------------------------------------
                    rhoUex[i+1][j] = rhoUex_out + rhoUex_in;
                    //std::cout << "j = " << j << " , "<<rhoUix_Rx << " , "<<rhoUex[i+1][j]<< " , "<<rhoUex_out << " , "<<rhoUex_in  << std::endl;
                    //------------------------------------

                    //wall-flux
                    //------------------------------------
                    rhoUex_wall[i+1][j] = 0.0;
                    //------------------------------------
                }
                //------------------------------------
            }
        }

        //for rhoUer
        //------------------------------------
        for (int iblock=0;iblock<2;iblock++){
            for (int i=i_flr_bl[iblock][0];i<=i_flr_bl[iblock][1];i++){ 
                for (int j=j_flr_bl[iblock][0];j<=j_flr_bl[iblock][1];j++){

                    //radial position at cell interface
                    //------------------------------------
                    double r_tmp = (r[j] + r[j-1])/2.0;
                    //------------------------------------

                    //moblity at cell interface
                    //------------------------------------
                    double mu_para_tmp = (mu_para[i][j] + mu_para[i][j-1])/2.0;
                    double mu_perp_tmp = (mu_perp[i][j] + mu_perp[i][j-1])/2.0;

                    double mu_star_tmp = mu_perp_tmp*c_mu;
                    //------------------------------------

                    //calculate inverse moblity tensor at cell interface
                    //------------------------------------
                    double mu_xx_tmp = (mu_xx[i][j] + mu_xx[i][j-1])/2.0;
                    double mu_xr_tmp = (mu_xr[i][j] + mu_xr[i][j-1])/2.0;
                    
                    double mu_inv_rr_tmp = mu_xx_tmp/(mu_para_tmp*mu_perp_tmp);
                    double mu_inv_xr_tmp = -mu_xr_tmp/(mu_para_tmp*mu_perp_tmp);
                    
                    //double mu_xx_tmp = (mu_xx[i][j] + mu_xx[i][j-1])/2.0;
                    //double mu_rr_tmp = (mu_rr[i][j] + mu_rr[i][j-1])/2.0;
                    //double mu_xr_tmp = (mu_xr[i][j] + mu_xr[i][j-1])/2.0;
                    //
                    //double deno = mu_xx_tmp*mu_rr_tmp - mu_xr_tmp*mu_xr_tmp;
                    //double mu_inv_rr_tmp = mu_xx_tmp/(deno + 1e-100);
                    //double mu_inv_xr_tmp = -mu_xr_tmp/(deno + 1e-100);
                    //------------------------------------
    
                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double Te_tmp = (Te[i][j] + Te[i][j-1])/2.0; //electron temperature
                    double dTedr_tmp = (Te[i][j] - Te[i][j-1])/dr*double(icon_dTe); //electron temperature gradient
                    double vd_tmp = mu_star_tmp*(-Er[i][j] - Boltz/e0*dTedr_tmp) + EPS; //electron drift velocity
                    double D_tmp = Boltz*Te_tmp/e0*mu_star_tmp; //electron diffusion coefficient
                    double z_tmp = vd_tmp*dr/D_tmp; //Peclet number
                    //------------------------------------

                    //interpolation of transverse flux component (should be old value)
                    //------------------------------------
                    double rhoUex_tmp = (r[j]*rhoUex_old_tmp[i+1][j] + r[j]*rhoUex_old_tmp[i][j]
                        + r[j-1]*rhoUex_old_tmp[i+1][j-1] + r[j-1]*rhoUex_old_tmp[i][j-1])/(4.0*r_tmp);
                    //------------------------------------

                    //update flux
                    //------------------------------------
                    rhoUer[i][j] =  
                          vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*rhoe[i][j-1]
                        - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*rhoe[i][j]
                        + rhoUer_old_tmp[i][j]
                        - mu_star_tmp*(mu_inv_xr_tmp*rhoUex_tmp + mu_inv_rr_tmp*rhoUer_old_tmp[i][j]);
                    //------------------------------------
                }
            }
        }
        //------------------------------------
        
        //lower-wall-BC (antenna)
        //------------------------------------
        for (int i=i_flc_bl[0][0];i<=i_flc_bl[0][1];i++){ 
            int j = j_flc_bl[0][0];

            //bulk-flux
            //------------------------------------
            rhoUer[i][j] = 0.0;
            //------------------------------------

            //wall-flux
            //------------------------------------
            double vth = sqrt(8.0*Boltz*Te[i][j]/(M_PI*masse));
            double delta_phi = fmax(-Er[i][j]*dr/2.0,0.0); //adjacent cell-center potential from wall
            double coef_angle = 1.0;
            double coef_mirror = 1.0;
            if(icon_eKineBC == 2) delta_phi = 0.0; //no need to extrapolate density
            if(icon_eFLRBC == 1){
                double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
                double angle = asin(fabs(Br[i][j]/Bmag)); //B-field angle from x-axis
                double vth2 = sqrt(2.0*Boltz*Te[i][j]/masse);
                double omega_ce = e0*Bmag/masse;
                double local_r = vth2/omega_ce;
                coef_angle = sqrt(pow(sin(angle),2) + pow(fmin(local_r/dx,1.0)*cos(angle),2));
            }else if(icon_eFLRBC == 2){
                double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
                double angle = asin(fabs(Br[i][j]/Bmag)); //B-field angle from x-axis
                double Halle = e0*Bmag/(masse*nu_m[i][j]);
                coef_angle = pow(sin(angle),2) + pow(cos(angle),2)/(1.0 + Halle*Halle);
            }
            if(icon_eMirBC == 1 && x[i]>=0.045 && x[i]<=0.055 && x[i]>=0.075){
                double lm = 0.015*M_PI;
                //double lm = M_PI/0.015*pow(x[i] - 0.065,2);
                double v_para = sqrt(2.0*Boltz*Te[i][j]/masse);
                double tau_m = lm/(2.0*v_para);
                coef_mirror = 1.0 - exp(-nu_m[i][j]*tau_m);
            }
            double GammaPerN_wall_Lr = 
                -0.25*vth*coef_angle*coef_mirror
                *exp(fmin(-e0*delta_phi/(Boltz*Te[i][j] + 1e-100),lim_exp))
                *(exp(-Uter[i][j]*Uter[i][j]) + uter[i][j]*sqrt(M_PI)*(erf(Uter[i][j]) - 1.0));

            rhoUer_wall[i][j] = rhoe[i][j]*GammaPerN_wall_Lr;
            //------------------------------------
        }
        //------------------------------------

        //centerline-normal-BC zero-flux
        //------------------------------------
        for (int i=i_flc_bl[1][0];i<=i_flc_bl[1][1];i++){ 
            int j = j_flc_bl[1][0];

            //bulk-flux
            //------------------------------------
            rhoUer[i][j] = 0.0;
            //------------------------------------

            //wall-flux
            //------------------------------------
            rhoUer_wall[i][j] = 0.0;
            //------------------------------------
        }
        //------------------------------------

        if(icon_open == 0){
            //upper-wall-flux
            //------------------------------------
            for (int i=i_flc_bl[0][0];i<=i_flc_bl[1][1];i++){ 
                int j = j_flc_bl[1][1];

                //bulk-flux
                //------------------------------------
                rhoUer[i][j+1] = 0.0;
                //------------------------------------

                //calulate kinetic vector flux
                //------------------------------------
                double vth = sqrt(8.0*Boltz*Te[i][j]/(M_PI*masse));
                double delta_phi =  fmax(Er[i][j+1]*dr/2.0,0.0); //adjacent cell-center potential from wall
                double coef_angle = 1.0;
                if(icon_eKineBC == 2) delta_phi = 0.0; //no need to extrapolate density
                if(icon_eFLRBC == 1){
                    double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
                    double angle = asin(fabs(Br[i][j]/Bmag)); //B-field angle from x-axis
                    double vth2 = sqrt(2.0*Boltz*Te[i][j]/masse);
                    double omega_ce = e0*Bmag/masse;
                    double local_r = vth2/omega_ce;
                    coef_angle = sqrt(pow(sin(angle),2) + pow(fmin(local_r/dx,1.0)*cos(angle),2));
                }else if(icon_eFLRBC == 2){
                    double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
                    double angle = asin(fabs(Br[i][j]/Bmag)); //B-field angle from x-axis
                    double Halle = e0*Bmag/(masse*nu_m[i][j]);
                    coef_angle = pow(sin(angle),2) + pow(cos(angle),2)/(1.0 + Halle*Halle);
                }
                double GammaPerN_wall_Rr = 
                    0.25*vth*coef_angle
                    *exp(fmin(-e0*delta_phi/(Boltz*Te[i][j] + 1e-100),lim_exp))
                    *(exp(-Uter[i][j+1]*Uter[i][j+1]) + uter[i][j+1]*sqrt(M_PI)*(erf(Uter[i][j+1]) + 1.0));

                rhoUer_wall[i][j+1] = rhoe[i][j]*GammaPerN_wall_Rr;
                //------------------------------------
            }
            //------------------------------------
        
        }else{
            if(icon_opGauss == 0){
                //upper-open-flux
                //------------------------------------
                for (int i=i_flc_bl[0][0];i<=i_flc_bl[1][1];i++){ 
                    int j = j_flc_bl[1][1];

                    //bulk-flux (equal to ion)
                    //------------------------------------
                    double rhoUir_Rr = rhoi[i][j]*Uir[i][j+1];
                    if(icon_ionConsv == 1){
                        rhoUir_Rr = rhoUir[i][j+1];
                    }
                    double Gamma_open_Rr = rhoUir_Rr;

                    rhoUer[i][j+1] =  Gamma_open_Rr;
                    //------------------------------------

                    //wall-flux
                    //------------------------------------
                    rhoUer_wall[i][j+1] = 0.0;
                    //------------------------------------
                }
                //------------------------------------
            }else{
                //upper-open-flux
                //------------------------------------
                for (int i=i_flc_bl[0][0];i<=i_flc_bl[1][1];i++){ 
                    int j = j_flc_bl[1][1];


                    //outgoing flux
                    //------------------------------------
                    double vth = sqrt(8.0*Boltz*Te[i][j]/(M_PI*masse));
                    double GammaPerN_wall_Rr = 
                        0.25*vth
                        *(exp(-uter[i][j+1]*uter[i][j+1]) + uter[i][j+1]*sqrt(M_PI)*(erf(uter[i][j+1]) + 1.0));

                    double rhoUer_out = rhoe[i][j]*GammaPerN_wall_Rr;
                    //------------------------------------

                    //ion bulk-flux
                    //------------------------------------
                    double rhoUir_Rr = rhoi[i][j]*Uir[i][j+1];
                    if(icon_ionConsv == 1){
                        rhoUir_Rr = rhoUir[i][j+1];
                    }
                    //------------------------------------

                    //incoming flux
                    //------------------------------------
                    double rhoUer_in = fmin(rhoUir_Rr - eps0/(e0*dt)*(1.0 - R_Eop)*Er_old[i][j+1] - rhoUer_out,0.0);
                    //double rhoUer_in = rhoUir_Rr - eps0/(e0*dt)*(1.0 - R_Eop)*Er_old[i][j+1] - rhoUer_out;
                    //------------------------------------

                    //bulk-flux
                    //------------------------------------
                    rhoUer[i][j+1] = rhoUer_out + rhoUer_in;
                    //------------------------------------

                    //wall-flux
                    //------------------------------------
                    rhoUer_wall[i][j+1] = 0.0;
                    //------------------------------------
                }
                //------------------------------------
            }
        }

        //Make coefficient for electron continuity equation for 5 point stencil
        //------------------------------------
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rR = (r[j]+r[j+1])/2.0;
                    double rL = (r[j]+r[j-1])/2.0;
                    //------------------------------------

                    //normalized radial position at upper and lower cell interface
                    //------------------------------------
                    double qR = (r[j]+r[j+1])/2.0/r[j];
                    double qL = (r[j]+r[j-1])/2.0/r[j];
                    //------------------------------------
                    
                    //moblity at cell interface
                    //------------------------------------
                    double mu_star_Lx = (mu_perp[i][j]+mu_perp[i-1][j])/2.0*c_mu;
                    double mu_star_Rx = (mu_perp[i][j]+mu_perp[i+1][j])/2.0*c_mu;
                    double mu_star_Lr = (mu_perp[i][j]+mu_perp[i][j-1])/2.0*c_mu;
                    double mu_star_Rr = (mu_perp[i][j]+mu_perp[i][j+1])/2.0*c_mu;
                    //------------------------------------

                    //electron temperature at cell interface
                    //------------------------------------
                    double Te_Lx = (Te[i][j]+Te[i-1][j])/2.0;
                    double Te_Rx = (Te[i][j]+Te[i+1][j])/2.0;
                    double Te_Lr = (Te[i][j]+Te[i][j-1])/2.0;
                    double Te_Rr = (Te[i][j]+Te[i][j+1])/2.0;
                    //------------------------------------

                    //electron diffusion coefficient at cell interface
                    //------------------------------------
                    double D_Lx = Boltz*Te_Lx/e0*mu_star_Lx;
                    double D_Rx = Boltz*Te_Rx/e0*mu_star_Rx;
                    double D_Lr = Boltz*Te_Lr/e0*mu_star_Lr;
                    double D_Rr = Boltz*Te_Rr/e0*mu_star_Rr;
                    //------------------------------------

                    //electron temperature gradient at cell interface
                    //------------------------------------
                    double dTedx_Lx = (Te[i  ][j  ] - Te[i-1][j  ])/dx*double(icon_dTe);
                    double dTedx_Rx = (Te[i+1][j  ] - Te[i  ][j  ])/dx*double(icon_dTe);
                    double dTedr_Lr = (Te[i  ][j  ] - Te[i  ][j-1])/dr*double(icon_dTe);
                    double dTedr_Rr = (Te[i  ][j+1] - Te[i  ][j  ])/dr*double(icon_dTe);
                    //------------------------------------

                    //electron drift velocity at cell interface
                    //------------------------------------
                    double vd_Lx = mu_star_Lx*(-Ex[i  ][j  ] - Boltz/e0*dTedx_Lx) + EPS;
                    double vd_Rx = mu_star_Rx*(-Ex[i+1][j  ] - Boltz/e0*dTedx_Rx) + EPS;
                    double vd_Lr = mu_star_Lr*(-Er[i  ][j  ] - Boltz/e0*dTedr_Lr) + EPS;
                    double vd_Rr = mu_star_Rr*(-Er[i  ][j+1] - Boltz/e0*dTedr_Rr) + EPS;
                    //------------------------------------

                    //Peclet number at cell interface
                    //------------------------------------
                    double z_Lx = vd_Lx*dx/D_Lx;
                    double z_Rx = vd_Rx*dx/D_Rx;
                    double z_Lr = vd_Lr*dr/D_Lr;
                    double z_Rr = vd_Rr*dr/D_Rr;
                    //------------------------------------

                    //BC - set adjacent flag (adjacent = 0, others = 1)
                    //adjacent "cell" to the wall for setting zero-flux on the wall
                    //------------------------------------
                    //left
                    double bLx_wall  = double(i!=i_flc_bl[0][0])*double(i!=i_flc_bl[1][0] || j>=j_flc_bl[0][0]);
                    double bLx  = bLx_wall;

                    //right
                    double bRx_open  = double(i!=i_flc_bl[1][1] || icon_open != 1);
                    double bRx_wall  = double(i!=i_flc_bl[1][1] || icon_open != 0);
                    double bRx  = bRx_open*bRx_wall;
                    
                    //lower
                    double bLr_wall  = double(j!=j_flc_bl[0][0] || i >=i_flc_bl[1][0]);
                    double bLr_cen = double(j!=j_flc_bl[1][0]);
                    double bLr  = bLr_wall*bLr_cen;

                    //upper
                    double bRr_open  = double(j!=j_flc_bl[1][1] || icon_open != 1);
                    double bRr_wall  = double(j!=j_flc_bl[1][1] || icon_open != 0);
                    double bRr  = bRr_open*bRr_wall;
                    //------------------------------------

                    //coefficient for electron continuity equation for 5 point stencil
                    //==============================================
                    //vertical and horizental 5 stencil
                    //------------------------------------
                    aP[i][j] =    vd_Rx/(1.0 - exp(fmin(-z_Rx,lim_exp)))*bRx*dr
                        +    vd_Lx/(exp(fmin( z_Lx,lim_exp)) - 1.0)*bLx*dr
                        + qR*vd_Rr/(1.0 - exp(fmin(-z_Rr,lim_exp)))*bRr*dx
                        + qL*vd_Lr/(exp(fmin( z_Lr,lim_exp)) - 1.0)*bLr*dx
                        + dx*dr/dt; //time-term
                    aE[i][j] =    vd_Rx/(exp(fmin( z_Rx,lim_exp)) - 1.0)*bRx*dr;
                    aW[i][j] =    vd_Lx/(1.0 - exp(fmin(-z_Lx,lim_exp)))*bLx*dr;
                    aN[i][j] = qR*vd_Rr/(exp(fmin( z_Rr,lim_exp)) - 1.0)*bRr*dx;
                    aS[i][j] = qL*vd_Lr/(1.0 - exp(fmin(-z_Lr,lim_exp)))*bLr*dx;
                    //------------------------------------

                    //RHS
                    //------------------------------------
                    double nabla_rhoUe_tmp = (rhoUex[i+1][j]*bRx_wall - rhoUex[i][j]*bLx_wall)/dx 
                                           + (qR*rhoUer[i][j+1]*bRr_wall - qL*rhoUer[i][j]*bLr_wall)/dr
                                           + (rhoUex_wall[i+1][j] - rhoUex_wall[i][j])/dx 
                                           + (qR*rhoUer_wall[i][j+1] - qL*rhoUer_wall[i][j])/dr;
                    b[i][j] = (-nabla_rhoUe_tmp 
                        + rate_ionize[i][j]
                        - (rhoe[i][j] - rhoe_old[i][j])/dt)*dx*dr; //time-term
                    //------------------------------------
                    //==============================================
                }
            }
        }
        //------------------------------------

        //Solver
        //------------------------------------
        if(icon_iter_rhoe == 0){
            //SOR
            //------------------------------------
            double alpha = 1.9;
            solver_SOR(aP,aE,aW,aN,aS,b,i_flc_bl,j_flc_bl,alpha,maxITR_SOR_rhoe,error_cnv_SOR_rhoe,0,drhoe);
            //------------------------------------
        }else{
            //Semi-coarse Multi-Grid by Hypre
            //------------------------------------
            solver_SMG(aP,aE,aW,aN,aS,b,i_flc_bl,j_flc_bl,maxITR_SOR_rhoe,error_cnv_SOR_rhoe,0,drhoe);
            //------------------------------------
        }
        //------------------------------------

        //Correction
        //==============================================
        //for rhoUex
        //------------------------------------
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flx_bl[iblock][0];i<=i_flx_bl[iblock][1];i++){
                for (int j=j_flx_bl[iblock][0];j<=j_flx_bl[iblock][1];j++){
                    
                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double mu_star_tmp = (mu_perp[i][j] + mu_perp[i-1][j])/2.0*c_mu;
                    double dTedx_tmp = (Te[i][j] - Te[i-1][j])/dx*double(icon_dTe);
                    double vd_tmp = mu_star_tmp*(-Ex[i][j] - Boltz/e0*dTedx_tmp) + EPS;
                    double Te_tmp = (Te[i][j] + Te[i-1][j])/2.0;
                    double D_tmp = Boltz*Te_tmp/e0*mu_star_tmp;
                    double z_tmp = vd_tmp*dx/D_tmp;
                    //------------------------------------

                    //calculate flux correction
                    //------------------------------------
                    double drhoUex = vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*drhoe[i-1][j]
                     - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*drhoe[i][j];
                    //------------------------------------
                    
                    //Update
                    //------------------------------------
                    rhoUex[i][j] = rhoUex[i][j] + drhoUex;
                    //------------------------------------
                }
            }
        }
        //------------------------------------

        //for rhoUer
        //------------------------------------
        for (int iblock=0;iblock<2;iblock++){
            for (int i=i_flr_bl[iblock][0];i<=i_flr_bl[iblock][1];i++){ 
                for (int j=j_flr_bl[iblock][0];j<=j_flr_bl[iblock][1];j++){
                    
                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double mu_star_tmp = (mu_perp[i][j] + mu_perp[i][j-1])/2.0*c_mu;
                    double dTedr_tmp = (Te[i  ][j  ] - Te[i  ][j-1])/dr*double(icon_dTe);
                    double vd_tmp = mu_star_tmp*(-Er[i][j] - Boltz/e0*dTedr_tmp) + EPS;
                    double Te_tmp = (Te[i][j] + Te[i][j-1])/2.0;
                    double D_tmp = Boltz*Te_tmp/e0*mu_star_tmp;
                    double z_tmp = vd_tmp*dr/D_tmp;
                    //------------------------------------

                    //calculate flux correction
                    //------------------------------------
                    double drhoUer =  vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*drhoe[i][j-1]
                     - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*drhoe[i][j];
                    //------------------------------------
                    
                    //Update
                    //------------------------------------
                    rhoUer[i][j] = rhoUer[i][j] + drhoUer;
                    //------------------------------------
                }
            }
        }
        //------------------------------------

        //for rhoe
        //------------------------------------
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                    //double rhoe_new = rhoe[i][j] + drhoe[i][j];
                    //if(rhoe_new > 0.0){
                    //    rhoe[i][j] = rhoe_new;
                    //}

                    rhoe[i][j] = rhoe[i][j] + drhoe[i][j];
                }
            }
        }
        //------------------------------------
        //==============================================

        //Caluculate residuals
        //=====================================================================================
        //rhoe
        //------------------------------------
        double error_rhoe2 = 0.0;
        if(icon_error2 == 0){ //normalized maximum
            for (int iblock=0;iblock<2;iblock++){ 
                for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                    for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                        double error_rhoe_tmp = fabs((rhoe[i][j]-rhoe_old_tmp[i][j])/(rhoe[i][j]+rhoe_old_tmp[i][j]+1e-100));
                        if(error_rhoe_tmp > error_rhoe2){
                            error_rhoe2 = error_rhoe_tmp;
                        }
                    }
                }
            }
        }else{ //normalized L2 norm
            double ncount_rhoe = 0;
            for (int iblock=0;iblock<2;iblock++){ 
                for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                    for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                        error_rhoe2 += pow(rhoe[i][j]-rhoe_old_tmp[i][j],2)/(pow(rhoe_old_tmp[i][j]+1e10,2));
                        ncount_rhoe++;
                    }
                }
            }
            error_rhoe2 = sqrt(error_rhoe2/double(ncount_rhoe));
        }
        //------------------------------------

        //x-directional flux
        //------------------------------------
        double error_rhoUex2 = 0.0;
        if(icon_error2 == 0){ //normalized maximum
            for (int iblock=0;iblock<2;iblock++){ 
                for (int i=i_flc_bl[iblock][0] + 1 - iblock;i<=i_flc_bl[iblock][1];i++){ //一番左右は更新しない
                    for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                        double error_rhoUex_tmp = fabs((rhoUex[i][j]-rhoUex_old_tmp[i][j])/(rhoUex[i][j]+rhoUex_old_tmp[i][j]+1e-100));
                        if(error_rhoUex_tmp > error_rhoUex2){
                            error_rhoUex2 = error_rhoUex_tmp;
                        }
                    }
                }
            }
        }else{ //normalized L2 norm
            double ncount_rhoUex = 0;
            for (int iblock=0;iblock<2;iblock++){ 
                for (int i=i_flc_bl[iblock][0] + 1 - iblock;i<=i_flc_bl[iblock][1];i++){ //一番左右は更新しない
                    for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                        error_rhoUex2 += pow(rhoUex[i][j]-rhoUex_old_tmp[i][j],2)/(pow(rhoUex_old_tmp[i][j]+1e16,2));
                        ncount_rhoUex++;
                    }
                }
            }
            error_rhoUex2 = sqrt(error_rhoUex2/double(ncount_rhoUex));
        }
        //------------------------------------

        //r-directional flux
        //------------------------------------
        double error_rhoUer2 = 0.0;
        if(icon_error2 == 0){ //normalized maximum
            for (int iblock=0;iblock<2;iblock++){
                for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                    for (int j=j_flc_bl[iblock][0] + 1;j<=j_flc_bl[iblock][1];j++){
                        double error_rhoUer_tmp = fabs((rhoUer[i][j]-rhoUer_old_tmp[i][j])/(rhoUer[i][j]+rhoUer_old_tmp[i][j]+1e-100));
                        if(error_rhoUer_tmp > error_rhoUer2){
                            error_rhoUer2 = error_rhoUer_tmp;
                        }
                    }
                }
            }
        }else{ //normalized L2 norm
            double ncount_rhoUer = 0;
            for (int iblock=0;iblock<2;iblock++){
                for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                    for (int j=j_flc_bl[iblock][0] + 1;j<=j_flc_bl[iblock][1];j++){
                        error_rhoUer2 += pow(rhoUer[i][j]-rhoUer_old_tmp[i][j],2)/(pow(rhoUer_old_tmp[i][j]+1e16,2))/CFL;
                        ncount_rhoUer++;
                    }
                }
            }
            error_rhoUer2 = sqrt(error_rhoUer2/double(ncount_rhoUer));
        }
        //------------------------------------
        //=====================================================================================
        
        double error_global = fmax(fmax(error_rhoe2,error_rhoUex2),error_rhoUer2);
        if(ncount % ndiv_out == 0){
            std::cout << "ncount = "<<ncount<<std::endl;
            std::cout << "error_rhoe = " << error_rhoe2
                <<  " error_rhoUex = " << error_rhoUex2
                <<  " error_rhoUer= " << error_rhoUer2
                <<  " error_max = " << fmax(fmax(error_rhoe2,error_rhoUex2),error_rhoUer2)
                << std::endl;
        }

        if(error_global < error_cnv_HES_rhoe || ncount >= maxITR_HES_rhoe){
            icon_end = 1;
        }

    } while (icon_end == 0);

    //Caluculate residuals
    //=====================================================================================
    //Electron density
    //------------------------------------
    error_rhoeps = 0.0;
    if(icon_error == 0){ //normalized maximum
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                    double error_tmp = fabs((rhoeps[i][j] - rhoeps_old[i][j])/(rhoeps[i][j] + rhoeps_old[i][j]+1e-100)*2.0/CFL);
                    if(error_tmp > error_rhoeps){
                        error_rhoeps = error_tmp;
                    }
                }
            }
        }
    }else{ //normalized L2 norm
        double ncount = 0;
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                    error_rhoeps += pow(rhoeps[i][j]-rhoeps_old[i][j],2)/(pow(rhoeps_old[i][j]+1e10,2))/CFL;
                    ncount++;
                }
            }
        }
        error_rhoeps = sqrt(error_rhoeps/double(ncount));
    }
    //------------------------------------

    //x-directional flux
    //------------------------------------
    error_Gx = 0.0;
    if(icon_error == 0){ //normalized maximum
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flx_bl[iblock][0];i<=i_flx_bl[iblock][1];i++){ 
                for (int j=j_flx_bl[iblock][0];j<=j_flx_bl[iblock][1];j++){
                    double error_tmp = fabs((Gx[i][j] - Gx_old[i][j])/(Gx[i][j] + Gx_old[i][j]+1e-100)*2.0/CFL);
                    if(error_tmp > error_Gx){
                        error_Gx = error_tmp;
                    }
                }
            }
        }
    }else{ //normalized L2 norm
        double ncount = 0;
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flx_bl[iblock][0];i<=i_flx_bl[iblock][1];i++){ 
                for (int j=j_flx_bl[iblock][0];j<=j_flx_bl[iblock][1];j++){
                    error_Gx += pow(Gx[i][j]-Gx_old[i][j],2)/(pow(Gx_old[i][j]+1.0,2))/CFL;
                    ncount++;
                }
            }
        }
        error_Gx = sqrt(error_Gx/double(ncount));
    }
    //------------------------------------

    //r-directional flux
    //------------------------------------
    error_Gr = 0.0;
    if(icon_error == 0){ //normalized maximum
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flr_bl[iblock][0];i<=i_flr_bl[iblock][1];i++){ 
                for (int j=j_flr_bl[iblock][0];j<=j_flr_bl[iblock][1];j++){
                    double error_tmp = fabs((Gr[i][j] - Gr_old[i][j])/(Gr[i][j] + Gr_old[i][j]+1e-100)*2.0/CFL);
                    if(error_tmp > error_Gr){
                        error_Gr = error_tmp;
                    }
                }
            }
        }
    }else{ //normalized L2 norm
        double ncount = 0;
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flr_bl[iblock][0];i<=i_flr_bl[iblock][1];i++){ 
                for (int j=j_flr_bl[iblock][0];j<=j_flr_bl[iblock][1];j++){
                    error_Gr += pow(Gr[i][j]-Gr_old[i][j],2)/(pow(Gr_old[i][j]+1.0,2)+1e-100)/CFL;
                    ncount++;
                }
            }
        }
        error_Gr = sqrt(error_Gr/double(ncount));
    }
    //------------------------------------
    //=====================================================================================

    //calculation of particle balance
    //------------------------------------
    if(itime%2000 == 0 || itime == ntime){
        
        //volume integrated values
        //------------------------------------
        double dNedt_sum         = 0.0;
        double nabla_rhoUe_sum   = 0.0;
        double rate_d_ionz_sum   = 0.0;
        double rate_s_ionz_sum   = 0.0;
        //------------------------------------
    
        //calculate volume integrated of each term
        //------------------------------------
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                    
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rR = (r[j]+r[j+1])/2.0;
                    double rL = (r[j]+r[j-1])/2.0;
                    //------------------------------------

                    //normalized radial position at upper and lower cell interface
                    //------------------------------------
                    double qR = (r[j]+r[j+1])/2.0/r[j];
                    double qL = (r[j]+r[j-1])/2.0/r[j];
                    //------------------------------------

                    //BC - set adjacent flag (adjacent = 0, others = 1)
                    //adjacent "cell" to the wall for setting zero-flux on the wall
                    //------------------------------------
                    //left
                    double bLx_wall  = double(i!=i_flc_bl[0][0])*double(i!=i_flc_bl[1][0] || j>=j_flc_bl[0][0]);
                    double bLx  = bLx_wall;

                    //right
                    double bRx_open  = double(i!=i_flc_bl[1][1] || icon_open != 1);
                    double bRx_wall  = double(i!=i_flc_bl[1][1] || icon_open != 0);
                    double bRx  = bRx_open*bRx_wall;
                    
                    //lower
                    double bLr_wall  = double(j!=j_flc_bl[0][0] || i >=i_flc_bl[1][0]);
                    double bLr_cen = double(j!=j_flc_bl[1][0]);
                    double bLr  = bLr_wall*bLr_cen;

                    //upper
                    double bRr_open  = double(j!=j_flc_bl[1][1] || icon_open != 1);
                    double bRr_wall  = double(j!=j_flc_bl[1][1] || icon_open != 0);
                    double bRr  = bRr_open*bRr_wall;
                    //------------------------------------

                    double nabla_rhoUe = (rhoUex[i+1][j]*bRx_wall - rhoUex[i][j]*bLx_wall)/dx 
                                       + (qR*rhoUer[i][j+1]*bRr_wall - qL*rhoUer[i][j]*bLr_wall)/dr
                                       + (rhoUex_wall[i+1][j] - rhoUex_wall[i][j])/dx 
                                       + (qR*rhoUer_wall[i][j+1] - qL*rhoUer_wall[i][j])/dr;
                    
                    double volume = 2.0*M_PI*r[j]*dr*dx;
                    
                    dNedt_sum         = dNedt_sum         + volume*(rhoe[i][j] - rhoe_old[i][j])/dt;
                    nabla_rhoUe_sum   = nabla_rhoUe_sum   + volume*nabla_rhoUe;
                    rate_d_ionz_sum   = rate_d_ionz_sum   + volume*rhoe[i][j]*nu_ionz[i][j];
                    rate_s_ionz_sum   = rate_s_ionz_sum   + volume*rhoe[i][j]*nu_ionzStep[i][j];
                }
            }
        }

        //output particle balance for electron
        //------------------------------------
        std::string char1="results/particle_balance_electron";
        std::string char2=std::to_string(nOut);
        std::string char_csv=".csv";
        std::ofstream outputfile8(char1+char2+char_csv);
        outputfile8 << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
        outputfile8<< "dN/dt (pcl/s), loss_to_bnd (pcl/s), direct-iz (pcl/s), sw-iz (pcl/s), sum (pcl/s)"<< std::endl;
        outputfile8<< dNedt_sum << "," << -nabla_rhoUe_sum << "," << rate_d_ionz_sum << ","<< rate_s_ionz_sum 
            << ","<< -dNedt_sum + (-nabla_rhoUe_sum + rate_d_ionz_sum + rate_s_ionz_sum)<< std::endl;
        outputfile8.close();
        //------------------------------------
    }
    //------------------------------------

    //for rhoUep
    //------------------------------------
    for (int iblock=0;iblock<2;iblock++){
        for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
            for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){

                //radial position at upper and lower cell interface
                //------------------------------------
                double rL = (r[j]+r[j-1])/2.0;
                double rR = (r[j]+r[j+1])/2.0;
                //------------------------------------

                double rhoUex_tmp = (   rhoUex[i+1][j] +    rhoUex[i][j])/2.0;
                double rhoUer_tmp = (rR*rhoUer[i][j+1] + rL*rhoUer[i][j])/(2.0*r[j]);

                if(icon_Bohm == 0){
                    rhoUep[i][j] = rhoUer_tmp*Omega_x[i][j] - rhoUex_tmp*Omega_r[i][j];
                }else{
                    rhoUep[i][j] = e0/(masse*nu_m[i][j])*(rhoUer_tmp*Bx[i][j] - rhoUex_tmp*Br[i][j]);
                }
            }
        }
    }

    //Don’t update unless the value is positive.
    //------------------------------------
    for (int iblock=0;iblock<2;iblock++){ 
        for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
            for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                if(rhoe[i][j]<0.0){
                    rhoe[i][j] = fmax(rhoe_old[i][j],0.0);
                }
            }
        }
    }
    //------------------------------------
    
    //Calculate kinetic energy for Te
    //------------------------------------
    for (int iblock=0;iblock<2;iblock++){
        for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
            for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                //radial position at upper and lower cell interface
                //------------------------------------
                double rL = (r[j]+r[j-1])/2.0;
                double rR = (r[j]+r[j+1])/2.0;
                //------------------------------------
                
                //normalized radial position at upper and lower cell interface
                //------------------------------------
                double qL = rL/r[j];
                double qR = rR/r[j];
                //------------------------------------

                //values for upwind discretization
                //*****************************************************
                //density at cell (rhoUer control volume) interface
                //------------------------------------
                double rhoUex_tmp = (   rhoUex[i+1][j] +    rhoUex[i][j])/2.0;
                double rhoUer_tmp = (rR*rhoUer[i][j+1] + rL*rhoUer[i][j])/(2.0*r[j]);
                //------------------------------------

                //velocity at cell (rhoUer control volume) interface
                //------------------------------------
                double Uex_tmp = rhoUex_tmp  /(rhoe[i][j] + rhoe_max*ratio_rhoe);
                double Uer_tmp = rhoUer_tmp  /(rhoe[i][j] + rhoe_max*ratio_rhoe);
                double Uep_tmp = rhoUep[i][j]/(rhoe[i][j] + rhoe_max*ratio_rhoe);
                //------------------------------------
       
                //kinetic energy density
                //------------------------------------
                rhokine[i][j] =  0.5*masse*(rhoUex_tmp*Uex_tmp + rhoUer_tmp*Uer_tmp + rhoUep[i][j]*Uep_tmp);
                //------------------------------------
            }
        }
    }

    
    //set values of ghost cells for ES-potential calculation
    //=====================================================================================
    
    //left-wall-BC (chamber)
    //------------------------------------
    for (int j=j_flc_bl[0][0];j<=j_flc_bl[0][1];j++){
        int i = i_flc_bl[0][0];
        rhoe[i-1][j] = rhoe[i][j];
    }
    //------------------------------------

    //left-wall-BC (magnet)
    //------------------------------------
    for (int j=j_flc_bl[1][0];j<=j_flc_bl[0][0]-1;j++){
        int i = i_flc_bl[1][0];
        rhoe[i-1][j] = rhoe[i][j];
    }
    //------------------------------------

    //lower-wall-BC (antenna)
    //------------------------------------
    for (int i=i_flc_bl[0][0];i<=i_flc_bl[0][1];i++){ 
        int j = j_flc_bl[0][0];
        rhoe[i][j-1] = rhoe[i][j];
    }

    //right-open-BC
    //------------------------------------
    for (int j=j_flc_bl[1][0];j<=j_flc_bl[1][1];j++){
        int i = i_flc_bl[1][1];
        rhoe[i+1][j] = rhoe[i][j];
    }
    //------------------------------------

    //upper-open-BC
    //------------------------------------
    for (int i=i_flc_bl[0][0];i<=i_flc_bl[1][1];i++){ 
        int j = j_flc_bl[0][1];
        rhoe[i][j+1] = rhoe[i][j];
    }
    //------------------------------------
    //=====================================================================================
    
    if(icon_open == 0 && icon_diele == 1){
        //calculate charge accumulation
        //------------------------------------
        for (int j=j_flc_bl[1][0];j<=j_flc_bl[1][1];j++){
            int i = i_flc_bl[1][1];

            double rR = (r[j]+r[j+1])/2.0;
            double rL = (r[j]+r[j-1])/2.0;
            double qR = rR/r[j];
            double qL = rL/r[j];

            double rhoi_Lx = (rhoi[i][j] + rhoi[i-1][j])/2.0;
            double rhoi_Rx = rhoi[i][j];
            double rhoi_Lr = (rhoi[i][j] + rhoi[i][j-1])/2.0;
            double rhoi_Rr = rhoi[i][j];

            scx[i+1][j] = scx_old[i+1][j] + dt*e0*(rhoUix_wall[i+1][j] - rhoUex_wall[i+1][j]);
            
            //std::cout << i << ","<< j << ","<<scx[i+1][j]
            //    << ","<< rhoUix_wall[i+1][j] << "," << rhoUex_wall[i+1][j] << std::endl;
        }

        for (int i=i_flc_bl[0][0];i<=i_flc_bl[1][1];i++){ 
            int j = j_flc_bl[0][1];
            double rR = (r[j]+r[j+1])/2.0;
            double rL = (r[j]+r[j-1])/2.0;
            double qR = rR/r[j];
            double qL = rL/r[j];

            double rhoi_Lx = (rhoi[i][j] + rhoi[i-1][j])/2.0;
            double rhoi_Rx = rhoi[i][j];
            double rhoi_Lr = (rhoi[i][j] + rhoi[i][j-1])/2.0;
            double rhoi_Rr = rhoi[i][j];

            scr[i][j+1] = scr_old[i][j+1] + dt*e0*(rhoUir_wall[i][j+1] - rhoUer_wall[i][j+1]);

            //std::cout << i << ","<< j << ","<<scr[i][j+1]
            //    << ","<< rhoUir_wall[i][j+1] << "," << rhoUer_wall[i][j+1] << std::endl;
        }
        //------------------------------------
    }
   

    //Caluculate residuals
    //=====================================================================================
    //Electron density
    //------------------------------------
    error_rhoe = 0.0;
    if(icon_error == 0){ //normalized maximum
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                    double error_tmp = fabs((rhoe[i][j] - rhoe_old[i][j])/(rhoe[i][j] + rhoe_old[i][j]+1e-100)*2.0/CFL);
                    if(error_tmp > error_rhoe){
                        error_rhoe = error_tmp;
                    }
                }
            }
        }
    }else{ //normalized L2 norm
        double ncount = 0;
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
                for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                    error_rhoe += pow(rhoe[i][j]-rhoe_old[i][j],2)/(pow(rhoe_old[i][j]+1e10,2))/CFL;
                    ncount++;
                }
            }
        }
        error_rhoe = sqrt(error_rhoe/double(ncount));
    }
    //------------------------------------

    //x-directional flux
    //------------------------------------
    error_rhoUex = 0.0;
    if(icon_error == 0){ //normalized maximum
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flx_bl[iblock][0];i<=i_flx_bl[iblock][1];i++){ 
                for (int j=j_flx_bl[iblock][0];j<=j_flx_bl[iblock][1];j++){
                    double error_tmp = fabs((rhoUex[i][j] - rhoUex_old[i][j])/(rhoUex[i][j] + rhoUex_old[i][j]+1e-100)*2.0/CFL);
                    if(error_tmp > error_rhoUex){
                        error_rhoUex = error_tmp;
                    }
                }
            }
        }
    }else{ //normalized L2 norm
        double ncount = 0;
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flx_bl[iblock][0];i<=i_flx_bl[iblock][1];i++){ 
                for (int j=j_flx_bl[iblock][0];j<=j_flx_bl[iblock][1];j++){
                    error_rhoUex += pow(rhoUex[i][j]-rhoUex_old[i][j],2)/(pow(rhoUex_old[i][j]+1e16,2))/CFL;
                    ncount++;
                }
            }
        }
        error_rhoUex = sqrt(error_rhoUex/double(ncount));
    }
    //------------------------------------

    //r-directional flux
    //------------------------------------
    error_rhoUer = 0.0;
    if(icon_error == 0){ //normalized maximum
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flr_bl[iblock][0];i<=i_flr_bl[iblock][1];i++){ 
                for (int j=j_flr_bl[iblock][0];j<=j_flr_bl[iblock][1];j++){
                    double error_tmp = fabs((rhoUer[i][j] - rhoUer_old[i][j])/(rhoUer[i][j] + rhoUer_old[i][j]+1e-100)*2.0/CFL);
                    if(error_tmp > error_rhoUer){
                        error_rhoUer = error_tmp;
                    }
                }
            }
        }
    }else{ //normalized L2 norm
        double ncount = 0;
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=i_flr_bl[iblock][0];i<=i_flr_bl[iblock][1];i++){ 
                for (int j=j_flr_bl[iblock][0];j<=j_flr_bl[iblock][1];j++){
                    error_rhoUer += pow(rhoUer[i][j]-rhoUer_old[i][j],2)/(pow(rhoUer_old[i][j]+1e16,2))/CFL;
                    ncount++;
                }
            }
        }
        error_rhoUer = sqrt(error_rhoUer/double(ncount));
    }
    //------------------------------------
    //=====================================================================================

}
