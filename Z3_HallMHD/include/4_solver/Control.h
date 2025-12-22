#pragma once
#include "0_basic/1_MPCNS_Parameter.h"
#include <math.h>
#include <iomanip>
class MHD_Control
{
public:
    /**
     * @brief   初始化，设置一些终止与控制的参数
     */
    void SetUp(Param *_par, int32_t number_of_var)
    {
        par = _par;
        residual.SetSize(number_of_var);
        residual_reference.SetSize(number_of_var);

        // 计算结束的判断标准
        Nstep = 0;
        maximum_Nstep = _par->GetInt("max_Nstep");
        tolerance = _par->GetDou("tolerance");
        maximum_Physic_Time = _par->GetDou("max_Time");
        // 计算中间控制参数
        Time_step = _par->GetDou("Time_step");
        CFL = _par->GetDou("CFL");
        Fileoutput_step = _par->GetInt("output_step");
        Resouput_step = _par->GetInt("output_residual");
        myid = _par->GetInt("myid");

        // 计时
        start = clock();
        start_program = start; // 记录开始时间
        cpu_time = 0.0;

        if (_par->GetBoo("continue_calc"))
        {
            Physic_Time = _par->GetDou("Physic_Time");
            Nstep = _par->GetInt("Nstep");

            residual_reference(0) = par->GetDou("Res_Ref_0");
            residual_reference(1) = par->GetDou("Res_Ref_1");
            residual_reference(2) = par->GetDou("Res_Ref_2");
            residual_reference(3) = par->GetDou("Res_Ref_3");
            residual_reference(4) = par->GetDou("Res_Ref_4");
            residual_reference(5) = par->GetDou("Res_Ref_5");
            residual_reference(6) = par->GetDou("Res_Ref_6");
            residual_reference(7) = par->GetDou("Res_Ref_7");
            residual = 1.0;
        }
        else
        {
            Physic_Time = 0.0;
            _par->AddParam("Physic_Time", Physic_Time);
        }
        _par->AddParam("Physic_Time_Step", 0.0);
    };

    /**
     * @brief   经过一个时间步后需要调用的更新control函数，记录步数，判断输出、终止程序
     * @remark  非定常计算才会存在物理时间，定常计算物理时间保持为0
     */
    void Update()
    {
        Nstep++;

        par->AddParam("Nstep", Nstep);
        par->AddParam("Physic_Time", Physic_Time);
        par->AddParam("Physic_Time_Step", Physic_Time_Step);

        //=========================================================================================
        // IF OUTPUT RESIDUAL
        //-----------------------------------------------
        residual_max = 0.0;
        for (int32_t i = 0; i < residual.Getsize1(); i++)
            residual_max = (residual_max < fabs(residual(i))) ? fabs(residual(i)) : residual_max;
        //-----------------------------------------------
        if_outres = (fmod(Nstep, Resouput_step) == 0 && myid == 0);
        if (if_outres)
        {
            //----------------------------------------------------------------------------
            end = clock();
            cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC; // 计算CPU执行时间，单位为秒
            start = end;
            //----------------------------------------------------------------------------
            std::cout << "Nstep=" << std::setw(10) << Nstep << "\tPhysicTime=" << std::setw(10) << std::setprecision(5) << Physic_Time
                      << "\t\tCPU_Time=" << std::setw(28) << std::setprecision(5) << cpu_time << std::endl;
            //   << "\tError=" << std::setw(10) << std::setprecision(5) << residual_max
            // << "\tCPU_Time=" << std::setw(10) << std::setprecision(5) << cpu_time << std::endl;
        }
        //=========================================================================================

        //=========================================================================================
        // IF OUTPUT FILE
        //-----------------------------------------------
        if_outfile = (fmod(Nstep, Fileoutput_step) == 0);
        //=========================================================================================

        //=========================================================================================
        // IF_STOP 判断是否终止程序计算过程，通过步数、残差以及最大时间同时判断
        //         非定常计算才会存在物理时间，定常计算物理时间保持为0
        //-----------------------------------------------
        if_stop = (Nstep >= maximum_Nstep);
        if_stop = if_stop || (Physic_Time >= maximum_Physic_Time);
        //-----------------------------------------------
        if_stop = if_stop || (residual_max < tolerance);
        //=========================================================================================
    };

public:
    Param *par;
    double1D residual, residual_reference;
    double Physic_Time, Physic_Time_Step, residual_max;

    // 计算中间控制参数
    double Time_step, CFL;
    int Fileoutput_step, Resouput_step, myid;

    // 终止程序的控制参数
    int Nstep, maximum_Nstep;
    double tolerance, maximum_Physic_Time;

    // 开关参数
    bool if_stop, if_outres, if_Calres, if_outfile;

    // 计时参数
    clock_t start, start_program, end;
    double cpu_time;

public:
    MHD_Control() {};
    ~MHD_Control() = default;
};
