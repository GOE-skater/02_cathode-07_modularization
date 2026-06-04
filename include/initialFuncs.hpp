#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>

#include "constants.hpp"
#include "params.hpp"
#include "miniFuncs.hpp"
#include "arrays.hpp"

using namespace std;

//*****************************************************************
//**                                                             **
//**           class InputModule                                 **
//**                                                             **
//*****************************************************************
class InitialFuncs
{
    private:

    public:
        void iniParam(Params &pm,GridCenter &gc,GridInterfaceX &gx,GridInterfaceR &gr);
        void makeBoundary(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);
};

//*****************************************************************
//**                                                             **
//**           void param                                        **
//**                                                             **
//*****************************************************************
void InitialFuncs::iniParam(Params &pm,GridCenter &gc,GridInterfaceX &gx,GridInterfaceR &gr)
{
    
    //******************* Initialization **************************
    //domain setting
    //--------------------------------
    pm.Lx = pm.xR - pm.xL;
    pm.Lr = pm.rmax - pm.rmin;
    pm.dx = pm.Lx/double(pm.ni);
    pm.dr = pm.Lr/double(pm.nj);
    std::cout << "dx = " << pm.dx<< " dr = " << pm.dr << std::endl;
    //--------------------------------

    //time step
    //--------------------------------
    pm.dt = pm.dt_ini;
    //--------------------------------

    //mesh Gegeration
    //--------------------------------
    for (int i = 0; i < pm.ni+2; i++){
        gc.x[i] = pm.xL+(i-0.5)*pm.dx;
    }
    for (int j = 0; j < pm.nj+2; j++){
        gc.r[j] = pm.rmin+(j-0.5)*pm.dr;
    }
    //--------------------------------

    //Give the r-directin current density Jr at x = 0.01 + dx/2, 0.003 < r < 0.006
    double i_exc = int(0.0002/pm.dx + 0.5); //5
    double j_exc_min = int(0.0008/pm.dr + 0.5) + 2; //6
    double j_exc_max = int(0.0020/pm.dr + 0.5); //10
    for (int j=j_exc_min;j<=j_exc_max;j++){
        //J1r[i_exc][j] = J1r_exc;

        double r_tmp = (gc.r[j] + gc.r[j-1])/2.0;
        gr.J1r[i_exc][j] = pm.J1r_exc/r_tmp/pm.dx;

        std::cout << i_exc << ","<<j << std::endl;
    }

    std::string folder_name = ("./results");
    mkdir(folder_name.c_str(), 0777);

}


//*****************************************************************
//**                                                             **
//**           void makeBoundary                                 **
//**                                                             **
//*****************************************************************
void InitialFuncs::makeBoundary(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr)
{
    
    double x_tmp = 0.0;
    double r_tmp = 0.0;
    
    //double x1 =  0.0e-3;
    //double x2 =  4.0e-3;
    //double x3 = 10.0e-3;
    //double x4 = 15.6e-3;
    //double x5 = 20.6e-3;
    //double x6 =  0.0e-3;

    pm.x1 =  0.0e-3;
    pm.x2 =  5.6e-3 + 7.8e-3; //5.5 mm
    pm.x3 = 10.0e-3 + 7.8e-3;
    pm.x4 = 15.0e-3 + 7.8e-3;
    pm.x5 = 20.0e-3 + 7.8e-3;
    pm.x6 =  0.0e-3 + 7.8e-3;
    
    pm.r1 =  0.8e-3;
    pm.r2 =  2.0e-3;
    pm.r3 =  4.0e-3;
    pm.r4 =  9.0e-3;
    pm.r5 =  2.0e-3;
    pm.r6 = 14.0e-3;

    //中性粒子の流入条件
    double Q_neutIn = pm.Q_neutIn_mgs*1e-6/pm.massi; //mg/s -> kg/s -> 個/s
    double S_in = 2.0*M_PI*pm.r4*pm.width_neutIn;
    pm.fn_In = Q_neutIn/S_in;
    
    gc.i_flc_bl[0][0] = int(pm.x6/pm.dx + 0.5) + 1; //1
    gc.j_flc_bl[0][0] = int(pm.r3/pm.dr + 0.5) + 1; //21
    gc.i_flc_bl[0][1] = int(pm.x2/pm.dx - 0.5) + 1; //20
    gc.j_flc_bl[0][1] = int(pm.r4/pm.dr - 0.5) + 1; //45

    gc.i_flc_bl[1][0] = int(pm.x2/pm.dx + 0.5) + 1; //21
    gc.j_flc_bl[1][0] = int(pm.r5/pm.dr + 0.5) + 1; //11
    gc.i_flc_bl[1][1] = int(pm.x4/pm.dx - 0.5) + 1; //75
    gc.j_flc_bl[1][1] = int(pm.r4/pm.dr - 0.5) + 1; //45

    gc.i_flc_bl[2][0] = int(pm.x2/pm.dx + 0.5) + 1; //21
    gc.j_flc_bl[2][0] = int(pm.r1/pm.dr + 0.5) + 1; //5
    gc.i_flc_bl[2][1] = int(pm.x3/pm.dx - 0.5) + 1; //50
    gc.j_flc_bl[2][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gc.i_flc_bl[3][0] = int(pm.x3/pm.dx + 0.5) + 1; //51
    gc.j_flc_bl[3][0] = int(0.0/pm.dr + 0.5) + 1; //1
    gc.i_flc_bl[3][1] = int(pm.x5/pm.dx - 0.5) + 1; //100
    gc.j_flc_bl[3][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gc.i_flc_bl[4][0] = int(pm.x5/pm.dx + 0.5) + 1; //101
    gc.j_flc_bl[4][0] = int(0.0/pm.dr + 0.5) + 1; //1
    gc.i_flc_bl[4][1] = int(pm.xR/pm.dx - 0.5) + 1; //158
    gc.j_flc_bl[4][1] = int(pm.rmax/pm.dr - 0.5) + 1; //70

    gc.i_flc_bl[5][0] = int(pm.x1/pm.dx + 0.5) + 1; //1
    gc.j_flc_bl[5][0] = int(pm.r1/pm.dr + 0.5) + 1; //5
    gc.i_flc_bl[5][1] = int(pm.x2/pm.dx - 0.5) + 1; //20
    gc.j_flc_bl[5][1] = int(pm.r2/pm.dr - 0.5) + 1; //10

    //std::cout << "flc" << std::endl;
    //std::cout << i_flc_bl[0][0] << " , "<<j_flc_bl[0][0] << " , "<< i_flc_bl[0][1] << " , "<<j_flc_bl[0][1] << std::endl;
    //std::cout << i_flc_bl[1][0] << " , "<<j_flc_bl[1][0] << " , "<< i_flc_bl[1][1] << " , "<<j_flc_bl[1][1] << std::endl;
    //std::cout << i_flc_bl[2][0] << " , "<<j_flc_bl[2][0] << " , "<< i_flc_bl[2][1] << " , "<<j_flc_bl[2][1] << std::endl;
    //std::cout << i_flc_bl[3][0] << " , "<<j_flc_bl[3][0] << " , "<< i_flc_bl[3][1] << " , "<<j_flc_bl[3][1] << std::endl;
    //std::cout << i_flc_bl[4][0] << " , "<<j_flc_bl[4][0] << " , "<< i_flc_bl[4][1] << " , "<<j_flc_bl[4][1] << std::endl;
    //std::cout << i_flc_bl[5][0] << " , "<<j_flc_bl[5][0] << " , "<< i_flc_bl[5][1] << " , "<<j_flc_bl[5][1] << std::endl;

    gx.i_flx_bl[0][0] = int(pm.x6/pm.dx + 0.5) + 2; //2
    gx.j_flx_bl[0][0] = int(pm.r3/pm.dr + 0.5) + 1; //21
    gx.i_flx_bl[0][1] = int(pm.x2/pm.dx - 0.5) + 2; //21
    gx.j_flx_bl[0][1] = int(pm.r4/pm.dr - 0.5) + 1; //47

    gx.i_flx_bl[1][0] = int(pm.x2/pm.dx + 0.5) + 2; //22
    gx.j_flx_bl[1][0] = int(pm.r5/pm.dr + 0.5) + 1; //11
    gx.i_flx_bl[1][1] = int(pm.x4/pm.dx - 0.5) + 1; //78
    gx.j_flx_bl[1][1] = int(pm.r4/pm.dr - 0.5) + 1; //47

    gx.i_flx_bl[2][0] = int(pm.x2/pm.dx + 0.5) + 2; //22
    gx.j_flx_bl[2][0] = int(pm.r1/pm.dr + 0.5) + 1; //5
    gx.i_flx_bl[2][1] = int(pm.x3/pm.dx - 0.5) + 2; //51
    gx.j_flx_bl[2][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gx.i_flx_bl[3][0] = int(pm.x3/pm.dx + 0.5) + 2; //52
    gx.j_flx_bl[3][0] = int(   0.0/pm.dr + 0.5) + 1; //1
    gx.i_flx_bl[3][1] = int(pm.x5/pm.dx - 0.5) + 2; //104
    gx.j_flx_bl[3][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gx.i_flx_bl[4][0] = int(pm.x5/pm.dx + 0.5) + 2; //105
    gx.j_flx_bl[4][0] = int(   0.0/pm.dr + 0.5) + 1; //1
    gx.i_flx_bl[4][1] = int(pm.xR/pm.dx - 0.5) + 1; //158
    gx.j_flx_bl[4][1] = int(pm.rmax/pm.dr - 0.5) + 1; //70

    gx.i_flx_bl[5][0] = int(pm.x1/pm.dx + 0.5) + 2; //2
    gx.j_flx_bl[5][0] = int(pm.r1/pm.dr + 0.5) + 1; //5
    gx.i_flx_bl[5][1] = int(pm.x2/pm.dx - 0.5) + 2; //21
    gx.j_flx_bl[5][1] = int(pm.r2/pm.dr - 0.5) + 1; //10

    //std::cout << "flx" << std::endl;
    //std::cout << i_flx_bl[0][0] << " , "<<j_flx_bl[0][0] << " , "<< i_flx_bl[0][1] << " , "<<j_flx_bl[0][1] << std::endl;
    //std::cout << i_flx_bl[1][0] << " , "<<j_flx_bl[1][0] << " , "<< i_flx_bl[1][1] << " , "<<j_flx_bl[1][1] << std::endl;
    //std::cout << i_flx_bl[2][0] << " , "<<j_flx_bl[2][0] << " , "<< i_flx_bl[2][1] << " , "<<j_flx_bl[2][1] << std::endl;
    //std::cout << i_flx_bl[3][0] << " , "<<j_flx_bl[3][0] << " , "<< i_flx_bl[3][1] << " , "<<j_flx_bl[3][1] << std::endl;
    //std::cout << i_flx_bl[4][0] << " , "<<j_flx_bl[4][0] << " , "<< i_flx_bl[4][1] << " , "<<j_flx_bl[4][1] << std::endl;
    //std::cout << i_flx_bl[5][0] << " , "<<j_flx_bl[5][0] << " , "<< i_flx_bl[5][1] << " , "<<j_flx_bl[5][1] << std::endl;

    gr.i_flr_bl[0][0] = int(pm.x6/pm.dx + 0.5) + 1; //1
    gr.j_flr_bl[0][0] = int(pm.r3/pm.dr + 0.5) + 2; //22
    gr.i_flr_bl[0][1] = int(pm.x2/pm.dx - 0.5) + 1; //20
    gr.j_flr_bl[0][1] = int(pm.r4/pm.dr - 0.5) + 1; //47

    gr.i_flr_bl[1][0] = int(pm.x2/pm.dx + 0.5) + 1; //21
    gr.j_flr_bl[1][0] = int(pm.r5/pm.dr + 0.5) + 1; //11
    gr.i_flr_bl[1][1] = int(pm.x4/pm.dx - 0.5) + 1; //78
    gr.j_flr_bl[1][1] = int(pm.r4/pm.dr - 0.5) + 1; //47

    gr.i_flr_bl[2][0] = int(pm.x2/pm.dx + 0.5) + 1; //21
    gr.j_flr_bl[2][0] = int(pm.r1/pm.dr + 0.5) + 2; //6
    gr.i_flr_bl[2][1] = int(pm.x3/pm.dx - 0.5) + 1; //50
    gr.j_flr_bl[2][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gr.i_flr_bl[3][0] = int(pm.x3/pm.dx + 0.5) + 1; //51
    gr.j_flr_bl[3][0] = int(   0.0/pm.dr + 0.5) + 2; //2
    gr.i_flr_bl[3][1] = int(pm.x5/pm.dx - 0.5) + 1; //103
    gr.j_flr_bl[3][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gr.i_flr_bl[4][0] = int(pm.x5/pm.dx + 0.5) + 1; //104
    gr.j_flr_bl[4][0] = int(   0.0/pm.dr + 0.5) + 2; //2
    gr.i_flr_bl[4][1] = int(pm.xR/pm.dx - 0.5) + 1; //158
    gr.j_flr_bl[4][1] = int(pm.rmax/pm.dr - 0.5) + 1; //70

    gr.i_flr_bl[5][0] = int(pm.x1/pm.dx + 0.5) + 1; //1
    gr.j_flr_bl[5][0] = int(pm.r1/pm.dr + 0.5) + 2; //6
    gr.i_flr_bl[5][1] = int(pm.x2/pm.dx - 0.5) + 1; //20
    gr.j_flr_bl[5][1] = int(pm.r2/pm.dr - 0.5) + 1; //10

    //std::cout << "flr" << std::endl;
    //std::cout << i_flr_bl[0][0] << " , "<<j_flr_bl[0][0] << " , "<< i_flr_bl[0][1] << " , "<<j_flr_bl[0][1] << std::endl;
    //std::cout << i_flr_bl[1][0] << " , "<<j_flr_bl[1][0] << " , "<< i_flr_bl[1][1] << " , "<<j_flr_bl[1][1] << std::endl;
    //std::cout << i_flr_bl[2][0] << " , "<<j_flr_bl[2][0] << " , "<< i_flr_bl[2][1] << " , "<<j_flr_bl[2][1] << std::endl;
    //std::cout << i_flr_bl[3][0] << " , "<<j_flr_bl[3][0] << " , "<< i_flr_bl[3][1] << " , "<<j_flr_bl[3][1] << std::endl;
    //std::cout << i_flr_bl[4][0] << " , "<<j_flr_bl[4][0] << " , "<< i_flr_bl[4][1] << " , "<<j_flr_bl[4][1] << std::endl;
    //std::cout << i_flr_bl[5][0] << " , "<<j_flr_bl[5][0] << " , "<< i_flr_bl[5][1] << " , "<<j_flr_bl[5][1] << std::endl;

    for (int iblock=0;iblock<5;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                gc.jdgBnd_flc[i][j] = 1;
            }
        }
    }

    for (int iblock=0;iblock<5;iblock++){ 
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){ 
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                gx.jdgBnd_flx[i][j] = 1;
            }
        }
    }

    for (int iblock=0;iblock<5;iblock++){ 
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                gr.jdgBnd_flr[i][j] = 1;
            }
        }
    }

    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                gc.jdgBnd_Ep[i][j] = 1;
            }
        }
    }

    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){ 
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                gx.jdgBnd_Ex[i][j] = 1;
            }
        }
    }

    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                gr.jdgBnd_Er[i][j] = 1;
            }
        }
    }

    //変換係数作成
    int kx_tmp = 0;
    int kr_tmp = 0;
    int kp_tmp = 0;
    int kfc_tmp = 0;
    int kfx_tmp = 0;
    int kfr_tmp = 0;

    for (int i=0;i<=ni+1;i++){
        for (int j=0;j<=nj+1;j++){
            if(jdgBnd_Ex[i][j]==1){
                kx[i][j] = kx_tmp;
                kx_tmp++;
                ikx.push_back(i);
                jkx.push_back(j);
            }

            if(jdgBnd_Er[i][j]==1){
                kr[i][j] = kr_tmp;
                kr_tmp++;
                ikr.push_back(i);
                jkr.push_back(j);
            }

            if(jdgBnd_Ep[i][j]==1){
                kp[i][j] = kp_tmp;
                kp_tmp++;
                ikp.push_back(i);
                jkp.push_back(j);
            }
        }
    }


    nkx = ikx.size();
    nkr = ikr.size();
    nkp = ikp.size();
    nk = nkx + nkr + nkp;
    std::cout << "nkx = " << nkx << " nkr = " << nkr << " nkp = " << nkp << " nk = " << nk  << std::endl;

    //Output boundary check file
    if(icon_chk == 1){

        std::ofstream outputfile1("results/boundary0.csv");
        outputfile1<<"i,j,x,r,jdgBnd_Ep,jdgBnd_Ex,jdgBnd_Er,jdgBnd_flc,jdgBnd_flx,jdgBnd_flr,kx,kr,kp,zero" << std::endl;

        for(int i=0;i<=ni+1;i++){
            for(int j=0;j<=nj+1;j++){
                outputfile1<< i << ","<< j << "," << x[i]<< ","<< r[j]
                    << "," << jdgBnd_Ep[i][j] << "," << jdgBnd_Ex[i][j]<< "," << jdgBnd_Er[i][j]
                    << "," << jdgBnd_flc[i][j]<< "," << jdgBnd_flx[i][j]<< "," << jdgBnd_flr[i][j]
                    << "," << kx[i][j]<< "," << kr[i][j]<< "," << kp[i][j]
                    << "," << 0.0<< std::endl;
            }
        }

        outputfile1.close();
    }
    
    //**************** 境界条件保持配列作成 ****************
    //********* Ex境界条件 (左) *********
    //z0
    for (int j=j_flc_bl[0][0]+1;j<=j_flc_bl[0][1]-1;j++){
        int i=i_flc_bl[0][0]+1;
        iBndWx.push_back(i);
        jBndWx.push_back(j);
        sBndWx.push_back(2); //ガウス
    }
    //z0
    {
        int i=i_flc_bl[0][0]+1;
        int j=j_flc_bl[0][0];
        iBndWx.push_back(i);
        jBndWx.push_back(j);
        sBndWx.push_back(3); //凹角 下側壁
    }
    //z0
    {
        int i=i_flc_bl[0][0]+1;
        int j=j_flc_bl[1][1];
        iBndWx.push_back(i);
        jBndWx.push_back(j);
        sBndWx.push_back(4); //凹角 上側壁
    }
    //z1
    for (int j=j_flc_bl[5][1]+1;j<=j_flc_bl[0][0]-1;j++){
        int i=i_flc_bl[2][0]+1;
        iBndWx.push_back(i);
        jBndWx.push_back(j);
        sBndWx.push_back(2); //ガウス
    }
    //z2
    for (int j=j_flc_bl[3][0];j<=j_flc_bl[2][0]-1;j++){
        int i=i_flc_bl[3][0]+1;
        iBndWx.push_back(i);
        jBndWx.push_back(j);
        sBndWx.push_back(2); //ガウス
    }
    //z4
    for (int j=j_flc_bl[1][0];j<=j_flc_bl[4][1]-1;j++){
        int i=i_flc_bl[4][0]+1;
        iBndWx.push_back(i);
        jBndWx.push_back(j);
        sBndWx.push_back(2); //ガウス
    }
    //z4
    {
        int i=i_flc_bl[4][0]+1;
        int j=j_flc_bl[4][1];
        iBndWx.push_back(i);
        jBndWx.push_back(j);
        sBndWx.push_back(4); //凹角 上側壁
    }
    //z6
    for (int j=j_flc_bl[5][0];j<=j_flc_bl[5][1];j++){
        int i=i_flc_bl[5][0]+1;
        iBndWx.push_back(i);
        jBndWx.push_back(j);
        sBndWx.push_back(0); //ディリクレ
    }
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileWx("results/BCx0.csv");
        outputfileWx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndWx.size();k++){
            int i = iBndWx[k];
            int j = jBndWx[k];
            
            outputfileWx << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndWx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileWx.close();
    }
    //********* Ex境界条件 (右) *********
    //z3
    for (int j=j_flc_bl[1][0];j<=j_flc_bl[1][1]-1;j++){
        int i=i_flc_bl[1][1];
        iBndEx.push_back(i);
        jBndEx.push_back(j);
        sBndEx.push_back(2); //ガウス
    }
    //z3
    {
        int i=i_flc_bl[1][1];
        int j=j_flc_bl[1][1];
        iBndEx.push_back(i);
        jBndEx.push_back(j);
        sBndEx.push_back(4); //凹角 上側壁
    }
    //z5
    for (int j=j_flc_bl[4][0];j<=j_flc_bl[4][1]-1;j++){
        int i=i_flc_bl[4][1];
        iBndEx.push_back(i);
        jBndEx.push_back(j);
        sBndEx.push_back(2); //ガウス
    }
    //z5
    {
        int i=i_flc_bl[4][1];
        int j=j_flc_bl[4][1];
        iBndEx.push_back(i);
        jBndEx.push_back(j);
        sBndEx.push_back(4); //凹角 上側壁
    }
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileEx("results/BCx1.csv");
        outputfileEx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndEx.size();k++){
            int i = iBndEx[k];
            int j = jBndEx[k];
            
            outputfileEx << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndEx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileEx.close();
    }
    //********* Ex境界条件 (下) *********
    //x0
    for (int i=i_flc_bl[0][0]+1;i<=i_flc_bl[0][1];i++){
        int j=j_flc_bl[0][0];
        iBndSx.push_back(i);
        jBndSx.push_back(j);
        sBndSx.push_back(0); //ディリクレ
    }
    //x0
    {
        int i=i_flc_bl[0][1]+1;
        int j=j_flc_bl[0][0];
        iBndSx.push_back(i);
        jBndSx.push_back(j);
        sBndSx.push_back(6); //凸角 右側Open
    }
    //x2
    for (int i=i_flc_bl[5][0]+1;i<=i_flc_bl[2][1];i++){
        int j=j_flc_bl[5][0];
        iBndSx.push_back(i);
        jBndSx.push_back(j);
        sBndSx.push_back(0); //ディリクレ
    }
    //x2
    {
        int i=i_flc_bl[2][1]+1;
        int j=j_flc_bl[5][0];
        iBndSx.push_back(i);
        jBndSx.push_back(j);
        sBndSx.push_back(6); //凸角 右側Open
    }
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileSx("results/BCx2.csv");
        outputfileSx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndSx.size();k++){
            int i = iBndSx[k];
            int j = jBndSx[k];
            
            outputfileSx << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndSx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileSx.close();
    }
    //********* Ex境界条件 (上) *********
    //x1
    for (int i=i_flc_bl[0][0]+1;i<=i_flc_bl[1][1];i++){
        int j=j_flc_bl[0][1];
        iBndNx.push_back(i);
        jBndNx.push_back(j);
        sBndNx.push_back(0); //ディリクレ
    }
    //x3
    for (int i=i_flc_bl[5][0]+1;i<=i_flc_bl[5][1];i++){
        int j=j_flc_bl[5][1];
        iBndNx.push_back(i);
        jBndNx.push_back(j);
        sBndNx.push_back(0); //ディリクレ
    }
    //x3
    {
        int i=i_flc_bl[5][1]+1;
        int j=j_flc_bl[5][1];
        iBndNx.push_back(i);
        jBndNx.push_back(j);
        sBndNx.push_back(6); //凸角 右側Open
    }
    //x4
    for (int i=i_flc_bl[1][1]+2;i<=i_flc_bl[3][1];i++){
        int j=j_flc_bl[3][1];
        iBndNx.push_back(i);
        jBndNx.push_back(j);
        sBndNx.push_back(0); //ディリクレ
    }
    //x4
    {
        int i=i_flc_bl[1][1]+1;
        int j=j_flc_bl[3][1];
        iBndNx.push_back(i);
        jBndNx.push_back(j);
        sBndNx.push_back(5); //凸角 左側Open
    }
    //x4
    {
        int i=i_flc_bl[3][1]+1;
        int j=j_flc_bl[3][1];
        iBndNx.push_back(i);
        jBndNx.push_back(j);
        sBndNx.push_back(6); //凸角 右側Open
    }
    //x5
    //for (int i=i_flc_bl[4][0]+1;i<=i_flc_bl[4][1];i++){
    //    int j=j_flc_bl[4][1];
    //    iBndNx.push_back(i);
    //    jBndNx.push_back(j);
    //    sBndNx.push_back(1); //開放
    //}
    //x5 (仮)
    for (int i=i_flc_bl[4][0]+1;i<=i_flc_bl[4][1];i++){
        int j=j_flc_bl[4][1];
        iBndNx.push_back(i);
        jBndNx.push_back(j);
        sBndNx.push_back(0); //ディリクレ
    }

    //出力
    if(icon_chk == 1){
        std::ofstream outputfileNx("results/BCx3.csv");
        outputfileNx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndNx.size();k++){
            int i = iBndNx[k];
            int j = jBndNx[k];
            
            outputfileNx << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndNx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileNx.close();
    }

    //********* Er境界条件 (左) *********
    //z0
    for (int j=j_flc_bl[0][0]+1;j<=j_flc_bl[0][1];j++){
        int i=i_flc_bl[0][0];
        iBndWr.push_back(i);
        jBndWr.push_back(j);
        sBndWr.push_back(0); //ディリクレ
    }
    //z1
    for (int j=j_flc_bl[5][1]+2;j<=j_flc_bl[0][0]-1;j++){
        int i=i_flc_bl[1][0];
        iBndWr.push_back(i);
        jBndWr.push_back(j);
        sBndWr.push_back(0); //ディリクレ
    }
    //z1
    {
        int i=i_flc_bl[1][0];
        int j=j_flc_bl[5][1]+1;
        iBndWr.push_back(i);
        jBndWr.push_back(j);
        sBndWr.push_back(5); //凸角 下側Open
    }
    //z1
    {
        int i=i_flc_bl[1][0];
        int j=j_flc_bl[0][0];
        iBndWr.push_back(i);
        jBndWr.push_back(j);
        sBndWr.push_back(6); //凸角 上側Open
    }
    //z2
    for (int j=j_flc_bl[3][0]+1;j<=j_flc_bl[2][0]-1;j++){
        int i=i_flc_bl[3][0];
        iBndWr.push_back(i);
        jBndWr.push_back(j);
        sBndWr.push_back(0); //ディリクレ
    }
    //z2
    {
        int i=i_flc_bl[3][0];
        int j=j_flc_bl[2][0];
        iBndWr.push_back(i);
        jBndWr.push_back(j);
        sBndWr.push_back(6); //凸角 上側Open
    }
    //z4
    for (int j=j_flc_bl[1][0]+1;j<=j_flc_bl[4][1];j++){
        int i=i_flc_bl[4][0];
        iBndWr.push_back(i);
        jBndWr.push_back(j);
        sBndWr.push_back(0); //ディリクレ
    }
    //z4
    {
        int i=i_flc_bl[4][0];
        int j=j_flc_bl[1][0];
        iBndWr.push_back(i);
        jBndWr.push_back(j);
        sBndWr.push_back(5); //凸角 下側Open
    }
    //z6
    for (int j=j_flc_bl[5][0]+1;j<=j_flc_bl[5][1];j++){
        int i=i_flc_bl[5][0];
        iBndWr.push_back(i);
        jBndWr.push_back(j);
        sBndWr.push_back(1); //開放
    }
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileWr("results/BCr0.csv");
        outputfileWr <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndWr.size();k++){
            int i = iBndWr[k];
            int j = jBndWr[k];
            
            outputfileWr << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndWr[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileWr.close();
    }
    //********* Er境界条件 (右) *********
    //z3
    for (int j=j_flc_bl[1][0]+1;j<=j_flc_bl[1][1];j++){
        int i=i_flc_bl[1][1];
        iBndEr.push_back(i);
        jBndEr.push_back(j);
        sBndEr.push_back(0); //ディリクレ
    }
    //z3
    {
        int i=i_flc_bl[1][1];
        int j=j_flc_bl[1][0];
        iBndEr.push_back(i);
        jBndEr.push_back(j);
        sBndEr.push_back(5); //凸角 下側Open
    }
    //z5
    for (int j=j_flc_bl[4][0]+1;j<=j_flc_bl[4][1];j++){
        int i=i_flc_bl[4][1];
        iBndEr.push_back(i);
        jBndEr.push_back(j);
        sBndEr.push_back(0); //ディリクレ
    }
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileEr("results/BCr1.csv");
        outputfileEr <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndEr.size();k++){
            int i = iBndEr[k];
            int j = jBndEr[k];
            
            outputfileEr << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndEr[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileEr.close();
    }
    //********* Er境界条件 (下) *********
    //x0
    for (int i=i_flc_bl[0][0]+1;i<=i_flc_bl[0][1];i++){
        int j=j_flc_bl[0][0]+1;
        iBndSr.push_back(i);
        jBndSr.push_back(j);
        sBndSr.push_back(2); //ガウス
    }
    //x0
    {
        int i=i_flc_bl[0][0];
        int j=j_flc_bl[0][0]+1;
        iBndSr.push_back(i);
        jBndSr.push_back(j);
        sBndSr.push_back(3); //凹角 左側壁
    }
    //x2
    for (int i=i_flc_bl[5][0];i<=i_flc_bl[2][1];i++){
        int j=j_flc_bl[5][0]+1;
        iBndSr.push_back(i);
        jBndSr.push_back(j);
        sBndSr.push_back(2); //ガウス
    }
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileSr("results/BCr2.csv");
        outputfileSr <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndSr.size();k++){
            int i = iBndSr[k];
            int j = jBndSr[k];
            
            outputfileSr << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndSr[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileSr.close();
    }
    //********* Er境界条件 (上) *********
    //x1
    for (int i=i_flc_bl[0][0]+1;i<=i_flc_bl[1][1]-1;i++){
        int j=j_flc_bl[0][1];
        iBndNr.push_back(i);
        jBndNr.push_back(j);
        sBndNr.push_back(2); //ガウス
    }
    //x1
    {
        int i=i_flc_bl[0][0];
        int j=j_flc_bl[0][1];
        iBndNr.push_back(i);
        jBndNr.push_back(j);
        sBndNr.push_back(3); //凹角 左側壁
    }
    //x1
    {
        int i=i_flc_bl[1][1];
        int j=j_flc_bl[0][1];
        iBndNr.push_back(i);
        jBndNr.push_back(j);
        sBndNr.push_back(4); //凹角 右側壁
    }
    //x3
    for (int i=i_flc_bl[5][0];i<=i_flc_bl[5][1];i++){
        int j=j_flc_bl[5][1];
        iBndNr.push_back(i);
        jBndNr.push_back(j);
        sBndNr.push_back(2); //ガウス
    }
    //x4
    for (int i=i_flc_bl[1][1]+1;i<=i_flc_bl[3][1];i++){
        int j=j_flc_bl[3][1];
        iBndNr.push_back(i);
        jBndNr.push_back(j);
        sBndNr.push_back(2); //ガウス
    }
    //x5
    for (int i=i_flc_bl[4][0]+1;i<=i_flc_bl[4][1]-1;i++){
        int j=j_flc_bl[4][1];
        iBndNr.push_back(i);
        jBndNr.push_back(j);
        sBndNr.push_back(2); //ガウス
    }
    //x5
    {
        int i=i_flc_bl[4][0];
        int j=j_flc_bl[4][1];
        iBndNr.push_back(i);
        jBndNr.push_back(j);
        sBndNr.push_back(3); //凹角 左側壁
    }
    //x5
    {
        int i=i_flc_bl[4][1];
        int j=j_flc_bl[4][1];
        iBndNr.push_back(i);
        jBndNr.push_back(j);
        sBndNr.push_back(4); //凹角 右側壁
    }
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileNr("results/BCr3.csv");
        outputfileNr <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndNr.size();k++){
            int i = iBndNr[k];
            int j = jBndNr[k];
            
            outputfileNr << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndNr[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileNr.close();
    }

    //********* Ep境界条件 (左) *********
    //z0
    for (int j=j_flc_bl[0][0];j<=j_flc_bl[0][1];j++){
        int i=i_flc_bl[0][0];
        iBndWp.push_back(i);
        jBndWp.push_back(j);
        sBndWp.push_back(0); //ディリクレ
    }
    //z1
    for (int j=j_flc_bl[5][1]+1;j<=j_flc_bl[0][0]-1;j++){
        int i=i_flc_bl[2][0];
        iBndWp.push_back(i);
        jBndWp.push_back(j);
        sBndWp.push_back(0); //ディリクレ
    }
    //z2
    for (int j=j_flc_bl[3][0];j<=j_flc_bl[2][0]-1;j++){
        int i=i_flc_bl[3][0];
        iBndWp.push_back(i);
        jBndWp.push_back(j);
        sBndWp.push_back(0); //ディリクレ
    }
    //z4
    for (int j=j_flc_bl[1][0];j<=j_flc_bl[4][1];j++){
        int i=i_flc_bl[4][0];
        iBndWp.push_back(i);
        jBndWp.push_back(j);
        sBndWp.push_back(0); //ディリクレ
    }
    //z6
    for (int j=j_flc_bl[5][0];j<=j_flc_bl[5][1];j++){
        int i=i_flc_bl[5][0];
        iBndWp.push_back(i);
        jBndWp.push_back(j);
        sBndWp.push_back(1); //開放
    }
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileWp("results/BCp0.csv");
        outputfileWp <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndWp.size();k++){
            int i = iBndWp[k];
            int j = jBndWp[k];
            
            outputfileWp << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndWp[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileWp.close();
    }
    //********* Ep境界条件 (右) *********
    //z3
    for (int j=j_flc_bl[1][0];j<=j_flc_bl[1][1];j++){
        int i=i_flc_bl[1][1];
        iBndEp.push_back(i);
        jBndEp.push_back(j);
        sBndEp.push_back(0); //ディリクレ
    }
    //z5
    for (int j=j_flc_bl[4][0];j<=j_flc_bl[4][1];j++){
        int i=i_flc_bl[4][1];
        iBndEp.push_back(i);
        jBndEp.push_back(j);
        sBndEp.push_back(0); //ディリクレ
    }
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileEp("results/BCp1.csv");
        outputfileEp <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndEp.size();k++){
            int i = iBndEp[k];
            int j = jBndEp[k];
            
            outputfileEp << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndEp[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileEp.close();
    }
    //********* Ep境界条件 (下) *********
    //x0
    for (int i=i_flc_bl[0][0];i<=i_flc_bl[0][1];i++){
        int j=j_flc_bl[0][0];
        iBndSp.push_back(i);
        jBndSp.push_back(j);
        sBndSp.push_back(0); //ディリクレ
    }
    //x2
    for (int i=i_flc_bl[5][0];i<=i_flc_bl[2][1];i++){
        int j=j_flc_bl[5][0];
        iBndSp.push_back(i);
        jBndSp.push_back(j);
        sBndSp.push_back(0); //ディリクレ
    }
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileSp("results/BCp2.csv");
        outputfileSp <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndSp.size();k++){
            int i = iBndSp[k];
            int j = jBndSp[k];
            
            outputfileSp << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndSp[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileSp.close();
    }
    //********* Ep境界条件 (上) *********
    //x1
    for (int i=i_flc_bl[0][0];i<=i_flc_bl[1][1];i++){
        int j=j_flc_bl[0][1];
        iBndNp.push_back(i);
        jBndNp.push_back(j);
        sBndNp.push_back(0); //ディリクレ
    }
    //x3
    for (int i=i_flc_bl[5][0];i<=i_flc_bl[5][1];i++){
        int j=j_flc_bl[5][1];
        iBndNp.push_back(i);
        jBndNp.push_back(j);
        sBndNp.push_back(0); //ディリクレ
    }
    //x4
    for (int i=i_flc_bl[1][1]+1;i<=i_flc_bl[3][1];i++){
        int j=j_flc_bl[3][1];
        iBndNp.push_back(i);
        jBndNp.push_back(j);
        sBndNp.push_back(0); //ディリクレ
    }
    //x5
    //for (int i=i_flc_bl[4][0];i<=i_flc_bl[4][1];i++){
    //    int j=j_flc_bl[4][1];
    //    iBndNp.push_back(i);
    //    jBndNp.push_back(j);
    //    sBndNp.push_back(1); //開放
    //}
    //x5 (仮)
    for (int i=i_flc_bl[4][0];i<=i_flc_bl[4][1];i++){
        int j=j_flc_bl[4][1];
        iBndNp.push_back(i);
        jBndNp.push_back(j);
        sBndNp.push_back(0); //ディリクレ
    }
    
    //出力
    if(icon_chk == 1){
        std::ofstream outputfileNp("results/BCp3.csv");
        outputfileNp <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<iBndNp.size();k++){
            int i = iBndNp[k];
            int j = jBndNp[k];
            
            outputfileNp << i << ","<< j 
                << ","<< x[i] << ","<<r[j] << ","<<sBndNp[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileNp.close();
    }
}
