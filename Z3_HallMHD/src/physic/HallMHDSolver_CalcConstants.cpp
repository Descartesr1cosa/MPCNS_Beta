#include "HallMHD_Solver.h"

void HallMHDSolver::calc_physical_constant(Param *par)
{
    List<double> temp;
    temp = par->GetDou_List("INITIAL");

    // ---- 常数 ----
    double gamma = par->GetDou_List("constant").data["gamma"];
    double NA = par->GetDou_List("constant").data["NA"];
    double R_uni = par->GetDou_List("constant").data["R_uni"];
    double q_e = par->GetDou_List("constant").data["q_e"];

    // ---- 流向单位向量 ----
    double c_y = temp.data["c_y"];
    double c_z = temp.data["c_z"];
    double c_x = std::sqrt(1.0 - c_y * c_y - c_z * c_z);

    // ---- 背景磁场（物理单位）----
    double Bx_phy = temp.data["Bx"];
    double By_phy = temp.data["By"];
    double Bz_phy = temp.data["Bz"];

    double B_ref = par->GetDou("B_ref");
    double L_ref = par->GetDou("L_ref");
    // 无量纲背景磁场
    double Bx = Bx_phy / B_ref;
    double By = By_phy / B_ref;
    double Bz = Bz_phy / B_ref;

    // 把无量纲 B_add 存进 Param，后面求解器/通量都可以直接用
    par->AddParam("B_add_x", Bx);
    par->AddParam("B_add_y", By);
    par->AddParam("B_add_z", Bz);

    // ---- 流体参考量 ----
    double U_ref = temp.data["U"];
    double n_ref = temp.data["n"];
    double T_ref = temp.data["T"];
    double Molecular_mass = temp.data["Molecular_mass"];

    double k_Boltz = R_uni / NA;
    double rho_ref = Molecular_mass / NA * n_ref;
    double p_ref = n_ref * k_Boltz * T_ref;

    double mu_mag = par->GetDou_List("constant").data["mu_mag"];
    double c_A = B_ref / std::sqrt(mu_mag * rho_ref);
    double Velocity_ref = U_ref;
    double M_A = Velocity_ref / c_A; // Alfven Mach 数
    double R_Larmor = Molecular_mass / NA * Velocity_ref / (B_ref * q_e);
    double phi = R_Larmor / L_ref;

    // 把无量纲 M_A^(-2)加入par便于调用
    par->AddParam("inver_MA2", 1.0 / (M_A * M_A));
    par->AddParam("phi", phi);

    //==============================================================================
    // 记录
    gamma_ = par_->GetDou_List("constant").data["gamma"];
    inver_MA2 = par_->GetDou("inver_MA2");
    hall_coef = inver_MA2 * par_->GetDou("phi");
}