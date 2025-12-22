#pragma once

#include <vector>

#include "0_basic/1_MPCNS_Parameter.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/2_MPCNS_Field.h"

#include "2_operators/CTOperators.h"
#include "util/Vector_tool.h"

class DipoleField
{
public:
    struct Dipole
    {
        Vec3 r0;  // 偶极子中心位置
        Vec3 m;   // 偶极矩向量（把常数系数吸收到m里）
        double a; // softening length（core）
    };

    // 偶极子磁场
    std::vector<Dipole> dips;

    void A_of_dipole(Vec3 &location, Vec3 &A)
    {
        A[0] = 0.0;
        A[1] = 0.0;
        A[2] = 0.0;
        for (auto &dip : dips)
        {
            Vec3 dr = location - dip.r0;
            double factor = dr.norm2() + dip.a * dip.a;
            factor = sqrt(factor) * factor;
            Vec3 A_temp = (dip.m ^ dr) / factor;
            A += A_temp;
        }
    }

    DipoleField() {};

    void load_from_param(Param *par)
    {
        dips.clear();

        const int n = par->GetInt("num_of_dipole");
        dips.reserve(n);

        double L_ref = par->GetDou("L_ref");
        double B_ref = par->GetDou("B_ref");
        double mu_mag = par->GetDou_List("constant").data["mu_mag"];

        double inver_m_mag_ref = 4 * 3.1415926535 * B_ref * L_ref * L_ref * L_ref / mu_mag;
        inver_m_mag_ref = 1.0 / inver_m_mag_ref; // 用于磁偶极子无量纲化的参数

        for (int i = 1; i <= n; ++i)
        {
            const std::string key = "Dipole" + std::to_string(i);

            auto dl = par->GetDou_List(key);

            const double x0 = dl.data["x0"];
            const double y0 = dl.data["y0"];
            const double z0 = dl.data["z0"];
            const double mx = dl.data["mx"];
            const double my = dl.data["my"];
            const double mz = dl.data["mz"];
            const double a = dl.data["a"];

            // --- 这里开始：把读到的数据映射到 Dipole 结构体 ---
            Dipole d;

            // 下面三行只是“占位映射示例”，你按实际含义改：
            d.r0 = Vec3(x0, y0, z0); // 有量纲长度
            d.m = Vec3(mx, my, mz);  // 有量纲磁偶极矩 A m m的量纲，后面无量纲化
            d.a = a;                 // 软化长度，有量纲

            d.r0 /= L_ref;
            d.a /= L_ref;
            d.m *= inver_m_mag_ref; // 无量纲化

            dips.push_back(d);
        }
    };

    void Build_Badd_FaceFlux(Grid *grd_, Field *fld_, Param *par_,
                             int Badd_xi_id, int Badd_eta_id, int Badd_zeta_id)
    {
        // 0) 清零
        //-----------------------------------------------------------------------------------------
        // Badd_ initialized as 0.0
        for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
        {
            // Badd_xi
            {
                auto &Badd_xi = fld_->field(Badd_xi_id, iblock);
                const Int3 &sub = Badd_xi.get_lo();
                const Int3 &sup = Badd_xi.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_xi(i, j, k, 0) = 0.0;
            }
            // Badd_eta
            {
                auto &Badd_eta = fld_->field(Badd_eta_id, iblock);
                const Int3 &sub = Badd_eta.get_lo();
                const Int3 &sup = Badd_eta.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_eta(i, j, k, 0) = 0.0;
            }
            // Badd_zeta
            {
                auto &Badd_zeta = fld_->field(Badd_zeta_id, iblock);
                const Int3 &sub = Badd_zeta.get_lo();
                const Int3 &sup = Badd_zeta.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_zeta(i, j, k, 0) = 0.0;
            }
        }

        // 1) 临时 edge 1-form 缓冲：A·dr（只在初始化函数里存在）
        std::vector<FieldBlock> Aadd_xi;
        std::vector<FieldBlock> Aadd_eta;
        std::vector<FieldBlock> Aadd_zeta;
        Aadd_xi.resize(fld_->num_blocks());
        Aadd_eta.resize(fld_->num_blocks());
        Aadd_zeta.resize(fld_->num_blocks());
        for (int b = 0; b < fld_->num_blocks(); ++b)
        {
            Aadd_xi[b].allocate(grd_->grids(b), {"Aadd_xi", StaggerLocation::EdgeXi, 1, 1});
            Aadd_eta[b].allocate(grd_->grids(b), {"Aadd_eta", StaggerLocation::EdgeEt, 1, 1});
            Aadd_zeta[b].allocate(grd_->grids(b), {"Aadd_zeta", StaggerLocation::EdgeZe, 1, 1});
        }

        for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
        {
            // Badd_xi
            {
                auto &A_xi = Aadd_xi[iblock];
                auto &x = grd_->grids(iblock).x;
                auto &y = grd_->grids(iblock).y;
                auto &z = grd_->grids(iblock).z;

                const Int3 &sub = A_xi.get_lo();
                const Int3 &sup = A_xi.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                        {
                            Vec3 x1 = {x(i, j, k), y(i, j, k), z(i, j, k)};
                            Vec3 x2 = {x(i + 1, j, k), y(i + 1, j, k), z(i + 1, j, k)};
                            Vec3 xm = 0.5 * (x1 + x2);
                            Vec3 dr = x2 - x1;
                            Vec3 A;
                            A_of_dipole(xm, A);        // A 是无量纲（因为 m 已无量纲化）
                            A_xi(i, j, k, 0) = A * dr; // A·dr，存到对应方向的 edge 标量里
                        }
            }
            // Badd_eta
            {
                auto &A_eta = Aadd_eta[iblock];
                auto &x = grd_->grids(iblock).x;
                auto &y = grd_->grids(iblock).y;
                auto &z = grd_->grids(iblock).z;

                const Int3 &sub = A_eta.get_lo();
                const Int3 &sup = A_eta.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                        {
                            Vec3 x1 = {x(i, j, k), y(i, j, k), z(i, j, k)};
                            Vec3 x2 = {x(i, j + 1, k), y(i, j + 1, k), z(i, j + 1, k)};
                            Vec3 xm = 0.5 * (x1 + x2);
                            Vec3 dr = x2 - x1;
                            Vec3 A;
                            A_of_dipole(xm, A);         // A 是无量纲（因为 m 已无量纲化）
                            A_eta(i, j, k, 0) = A * dr; // A·dr，存到对应方向的 edge 标量里
                        }
            }
            // Badd_zeta
            {
                auto &A_zeta = Aadd_zeta[iblock];
                auto &x = grd_->grids(iblock).x;
                auto &y = grd_->grids(iblock).y;
                auto &z = grd_->grids(iblock).z;

                const Int3 &sub = A_zeta.get_lo();
                const Int3 &sup = A_zeta.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                        {
                            Vec3 x1 = {x(i, j, k), y(i, j, k), z(i, j, k)};
                            Vec3 x2 = {x(i, j, k + 1), y(i, j, k + 1), z(i, j, k + 1)};
                            Vec3 xm = 0.5 * (x1 + x2);
                            Vec3 dr = x2 - x1;
                            Vec3 A;
                            A_of_dipole(xm, A);          // A 是无量纲（因为 m 已无量纲化）
                            A_zeta(i, j, k, 0) = A * dr; // A·dr，存到对应方向的 edge 标量里
                        }
            }
        }

        // 2) curl(edge)->face 得到偶极子部分的 Badd_face，然后累加到 IMF 上
        for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
        {
            auto &Badd_xi = fld_->field(Badd_xi_id, iblk);
            auto &Badd_eta = fld_->field(Badd_eta_id, iblk);
            auto &Badd_zeta = fld_->field(Badd_zeta_id, iblk);

            auto &A_xi = Aadd_xi[iblk];
            auto &A_eta = Aadd_eta[iblk];
            auto &A_zeta = Aadd_zeta[iblk];
            CTOperators::CurlEdgeToFace(iblk, A_xi, A_eta, A_zeta, Badd_xi, Badd_eta, Badd_zeta, /*multiper=*/1.0); //  B = curl A
        }

        // 3) IMF 直接加到 face 磁通：Phi = B0 · S_face
        //    这里 par 里有 Bx By Bz（无量纲或已归一）
        Vec3 B0(par_->GetDou("B_add_x"), par_->GetDou("B_add_y"), par_->GetDou("B_add_z"));
        for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
        {
            // Badd_xi
            {
                auto &Badd_xi = fld_->field(Badd_xi_id, iblock);
                auto &Aera = fld_->field("JDxi", iblock);
                const Int3 &sub = Badd_xi.get_lo();
                const Int3 &sup = Badd_xi.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_xi(i, j, k, 0) += (B0[0] * Aera(i, j, k, 0) + B0[1] * Aera(i, j, k, 1) + B0[2] * Aera(i, j, k, 2));
            }
            // Badd_eta
            {
                auto &Badd_eta = fld_->field(Badd_eta_id, iblock);
                auto &Aera = fld_->field("JDet", iblock);
                const Int3 &sub = Badd_eta.get_lo();
                const Int3 &sup = Badd_eta.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_eta(i, j, k, 0) += (B0[0] * Aera(i, j, k, 0) + B0[1] * Aera(i, j, k, 1) + B0[2] * Aera(i, j, k, 2));
            }
            // Badd_zeta
            {
                auto &Badd_zeta = fld_->field(Badd_zeta_id, iblock);
                auto &Aera = fld_->field("JDze", iblock);

                const Int3 &sub = Badd_zeta.get_lo();
                const Int3 &sup = Badd_zeta.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_zeta(i, j, k, 0) += (B0[0] * Aera(i, j, k, 0) + B0[1] * Aera(i, j, k, 1) + B0[2] * Aera(i, j, k, 2));
            }
        }
    };
};