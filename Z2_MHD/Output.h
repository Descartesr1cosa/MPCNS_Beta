#pragma once
#include "0_basic/DEFINE.h"
#include "0_basic/1_MPCNS_Parameter.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/2_MPCNS_Field.h"

class MHD_Output
{

public:
    std::string filename;
    Param *par;
    Field *fld;
    std::vector<std::string> var_defaut_plt_name;
    std::vector<std::string> var_defaut_bin_name;
    int32_t n_var_plt; // 变量个数，ZoneHeader会用

    MHD_Output() {};
    ~MHD_Output() = default;

    void SetUp(Param *_par, Field *_fld)
    {
        par = _par;
        fld = _fld;
        get_filename(par);
        get_iffile(fld);
        var_defaut_plt_name = {"rho", "u", "v", "w", "p", "Bx", "By", "Bz"};
        var_defaut_bin_name = {"U_", "B_xi", "B_eta", "B_zeta"};
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
    void output_bin_field() { output_bin_field(var_defaut_bin_name); };
    // 根据物理场名称列表输出
    void output_bin_field(const std::vector<std::string> &field_names);

    // 输出tecplot二进制文件，用于快速获取流场
    void output_plt_field() { output_plt_field(var_defaut_plt_name); };
    // 可以指定要输出的字段名
    void output_plt_field(const std::vector<std::string> &var_list);
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