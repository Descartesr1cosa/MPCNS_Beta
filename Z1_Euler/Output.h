#pragma once
#include "0_basic/DEFINE.h"
#include "0_basic/1_MPCNS_Parameter.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/2_MPCNS_Field.h"

class Euler_Output
{

public:
    std::string filename;
    Param *par;
    Field *fld;

    Euler_Output() {};
    ~Euler_Output() = default;

    void SetUp(Param *_par, Field *_fld)
    {
        par = _par;
        fld = _fld;
        get_filename(par);
        get_iffile(fld);
    }

    //=====================================================
    // 输出
    void output_field()
    {
        output_bin_field();
        output_plt_field();
    };

    //=====================================================
    // 输出二进制文件，用于续算需求
    void output_bin_field();
    // 输出tecplot二进制文件，用于快速获取流场
    void output_plt_field();
    //=====================================================

private:
    //=====================================================
    // Tecplot输出必备函数
    void plt_write_Headersection(const std::string &fieldname, std::vector<std::string> &var_name, FILE *file);
    void plt_write_str(const char *str, FILE *file);
    void plt_write_ZoneHeadersection(const std::string &Zone_name, int *mxyz, FILE *file);
    void plt_write_DataSection(int &n_var, FILE *file);
    //=====================================================
    // 输出ASCII格式的网格
    // void output_grids(Grid *grd, Param *par);
    // 获取本进程的文件名
    void get_filename(Param *par);
    // 获取本进程是否输出（多用于不同物理块耦合求解）
    void get_iffile(Field *fld) {};
    //=====================================================
};