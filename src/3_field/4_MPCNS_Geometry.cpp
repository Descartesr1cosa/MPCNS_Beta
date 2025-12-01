#include "3_field/2_MPCNS_Field.h"
#include "array"

void Field::build_geometry()
{
    register_field(FieldDescriptor{"Jac", StaggerLocation::Cell, 1, 0});
    register_field(FieldDescriptor{"JDxi", StaggerLocation::FaceXi, 3, 0});
    register_field(FieldDescriptor{"JDet", StaggerLocation::FaceEt, 3, 0});
    register_field(FieldDescriptor{"JDze", StaggerLocation::FaceZe, 3, 0});

    auto dot = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    auto cross = [&](const std::array<double, 3> &a, const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]};
    };
    auto minus = [&](const std::array<double, 3> &a, const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    };
    auto plus = [&](const std::array<double, 3> &a, const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    };

    {
        // Calc_ JDxi  = Jac\nabla\xi = Area
        auto &data_ = field("JDxi");
        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &data = data_[ib];
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            Int3 inner_range_lo = data.inner_lo();
            Int3 inner_range_hi = data.inner_hi();
            std::array<double, 3> ori, ori_xp, ori_yp, ori_xyp, rx, ry, rx_, ry_, Area;
            for (int i = inner_range_lo.i; i < inner_range_hi.i; i++)
                for (int j = inner_range_lo.j; j < inner_range_hi.j; j++)
                    for (int k = inner_range_lo.k; k < inner_range_hi.k; k++)
                    {
                        ori = {x(i, j, k), y(i, j, k), z(i, j, k)};
                        ori_xp = {x(i, j + 1, k), y(i, j + 1, k), z(i, j + 1, k)};
                        ori_yp = {x(i, j, k + 1), y(i, j, k + 1), z(i, j, k + 1)};
                        ori_xyp = {x(i, j + 1, k + 1), y(i, j + 1, k + 1), z(i, j + 1, k + 1)};
                        rx = minus(ori_xp, ori);
                        ry = minus(ori_yp, ori);
                        rx_ = minus(ori_yp, ori_xyp);
                        ry_ = minus(ori_xp, ori_xyp);
                        Area = plus(cross(rx, ry), cross(rx_, ry_));
                        data(i, j, k, 0) = 0.5 * Area[0];
                        data(i, j, k, 1) = 0.5 * Area[1];
                        data(i, j, k, 2) = 0.5 * Area[2];
                    }
        }
    }

    {
        // Calc_ JDet  = Jac\nabla\eta = Area
        auto &data_ = field("JDet");
        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &data = data_[ib];
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            Int3 inner_range_lo = data.inner_lo();
            Int3 inner_range_hi = data.inner_hi();
            std::array<double, 3> ori, ori_xp, ori_yp, ori_xyp, rx, ry, rx_, ry_, Area;
            for (int i = inner_range_lo.i; i < inner_range_hi.i; i++)
                for (int j = inner_range_lo.j; j < inner_range_hi.j; j++)
                    for (int k = inner_range_lo.k; k < inner_range_hi.k; k++)
                    {
                        ori = {x(i, j, k), y(i, j, k), z(i, j, k)};
                        ori_xp = {x(i, j, k + 1), y(i, j, k + 1), z(i, j, k + 1)};
                        ori_yp = {x(i + 1, j, k), y(i + 1, j, k), z(i + 1, j, k)};
                        ori_xyp = {x(i + 1, j, k + 1), y(i + 1, j, k + 1), z(i + 1, j, k + 1)};
                        rx = minus(ori_xp, ori);
                        ry = minus(ori_yp, ori);
                        rx_ = minus(ori_yp, ori_xyp);
                        ry_ = minus(ori_xp, ori_xyp);
                        Area = plus(cross(rx, ry), cross(rx_, ry_));
                        data(i, j, k, 0) = 0.5 * Area[0];
                        data(i, j, k, 1) = 0.5 * Area[1];
                        data(i, j, k, 2) = 0.5 * Area[2];
                    }
        }
    }

    {
        // Calc_ JDze  = Jac\nabla\zeta = Area
        auto &data_ = field("JDze");
        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &data = data_[ib];
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            Int3 inner_range_lo = data.inner_lo();
            Int3 inner_range_hi = data.inner_hi();
            std::array<double, 3> ori, ori_xp, ori_yp, ori_xyp, rx, ry, rx_, ry_, Area;
            for (int i = inner_range_lo.i; i < inner_range_hi.i; i++)
                for (int j = inner_range_lo.j; j < inner_range_hi.j; j++)
                    for (int k = inner_range_lo.k; k < inner_range_hi.k; k++)
                    {
                        ori = {x(i, j, k), y(i, j, k), z(i, j, k)};
                        ori_xp = {x(i + 1, j, k), y(i + 1, j, k), z(i + 1, j, k)};
                        ori_yp = {x(i, j + 1, k), y(i, j + 1, k), z(i, j + 1, k)};
                        ori_xyp = {x(i + 1, j + 1, k), y(i + 1, j + 1, k), z(i + 1, j + 1, k)};
                        rx = minus(ori_xp, ori);
                        ry = minus(ori_yp, ori);
                        rx_ = minus(ori_yp, ori_xyp);
                        ry_ = minus(ori_xp, ori_xyp);
                        Area = plus(cross(rx, ry), cross(rx_, ry_));
                        data(i, j, k, 0) = 0.5 * Area[0];
                        data(i, j, k, 1) = 0.5 * Area[1];
                        data(i, j, k, 2) = 0.5 * Area[2];
                    }
        }
    }

    {
        // Calc_ Jac = 1 / 3* \SUM Area \cdot dr
        auto &Jac_ = field("Jac");
        auto &Axi_ = field("JDxi");
        auto &Aeta_ = field("JDet");
        auto &Azeta_ = field("JDze");

        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &Jac = Jac_[ib];
            auto &A_xi = Axi_[ib];
            auto &A_et = Aeta_[ib];
            auto &A_ze = Azeta_[ib];

            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            auto &cx = grd->grids(ib).dual_x; // cell center
            auto &cy = grd->grids(ib).dual_y;
            auto &cz = grd->grids(ib).dual_z;

            Int3 lo = Jac.inner_lo();
            Int3 hi = Jac.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        std::array<double, 3> xc = {cx(i + 1, j + 1, k + 1),
                                                    cy(i + 1, j + 1, k + 1),
                                                    cz(i + 1, j + 1, k + 1)};

                        double V = 0.0;

                        auto dot = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
                        {
                            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                        };
                        auto minus = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
                        {
                            return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
                        };

                        // 下面要小心：对 minus 面，把 A 反号变成“对 cell 的外法向”
                        // 并且 face center 要用对应那 4 个顶点的平均

                        // --- xi- face at i ---
                        std::array<double, 3> Af_xm = {-A_xi(i, j, k, 0),
                                                       -A_xi(i, j, k, 1),
                                                       -A_xi(i, j, k, 2)};
                        std::array<double, 3> xf_xm = {
                            0.25 * (x(i, j, k) + x(i, j + 1, k) + x(i, j, k + 1) + x(i, j + 1, k + 1)),
                            0.25 * (y(i, j, k) + y(i, j + 1, k) + y(i, j, k + 1) + y(i, j + 1, k + 1)),
                            0.25 * (z(i, j, k) + z(i, j + 1, k) + z(i, j, k + 1) + z(i, j + 1, k + 1))};
                        V += dot(Af_xm, minus(xf_xm, xc));

                        // --- xi+ face at i+1 ---
                        std::array<double, 3> Af_xp = {A_xi(i + 1, j, k, 0),
                                                       A_xi(i + 1, j, k, 1),
                                                       A_xi(i + 1, j, k, 2)};
                        std::array<double, 3> xf_xp = {
                            0.25 * (x(i + 1, j, k) + x(i + 1, j + 1, k) + x(i + 1, j, k + 1) + x(i + 1, j + 1, k + 1)),
                            0.25 * (y(i + 1, j, k) + y(i + 1, j + 1, k) + y(i + 1, j, k + 1) + y(i + 1, j + 1, k + 1)),
                            0.25 * (z(i + 1, j, k) + z(i + 1, j + 1, k) + z(i + 1, j, k + 1) + z(i + 1, j + 1, k + 1))};
                        V += dot(Af_xp, minus(xf_xp, xc));

                        // --- eta- face at j ---
                        std::array<double, 3> Af_em = {-A_et(i, j, k, 0),
                                                       -A_et(i, j, k, 1),
                                                       -A_et(i, j, k, 2)};
                        std::array<double, 3> xf_em = {
                            0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j, k + 1) + x(i + 1, j, k + 1)),
                            0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j, k + 1) + y(i + 1, j, k + 1)),
                            0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j, k + 1) + z(i + 1, j, k + 1))};
                        V += dot(Af_em, minus(xf_em, xc));

                        // --- eta+ face at j+1 ---
                        std::array<double, 3> Af_ep = {A_et(i, j + 1, k, 0),
                                                       A_et(i, j + 1, k, 1),
                                                       A_et(i, j + 1, k, 2)};
                        std::array<double, 3> xf_ep = {
                            0.25 * (x(i, j + 1, k) + x(i + 1, j + 1, k) + x(i, j + 1, k + 1) + x(i + 1, j + 1, k + 1)),
                            0.25 * (y(i, j + 1, k) + y(i + 1, j + 1, k) + y(i, j + 1, k + 1) + y(i + 1, j + 1, k + 1)),
                            0.25 * (z(i, j + 1, k) + z(i + 1, j + 1, k) + z(i, j + 1, k + 1) + z(i + 1, j + 1, k + 1))};
                        V += dot(Af_ep, minus(xf_ep, xc));

                        // --- zeta- face at k ---
                        std::array<double, 3> Af_zm = {-A_ze(i, j, k, 0),
                                                       -A_ze(i, j, k, 1),
                                                       -A_ze(i, j, k, 2)};
                        std::array<double, 3> xf_zm = {
                            0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j + 1, k) + x(i + 1, j + 1, k)),
                            0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j + 1, k) + y(i + 1, j + 1, k)),
                            0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j + 1, k) + z(i + 1, j + 1, k))};
                        V += dot(Af_zm, minus(xf_zm, xc));

                        // --- zeta+ face at k+1 ---
                        std::array<double, 3> Af_zp = {A_ze(i, j, k + 1, 0),
                                                       A_ze(i, j, k + 1, 1),
                                                       A_ze(i, j, k + 1, 2)};
                        std::array<double, 3> xf_zp = {
                            0.25 * (x(i, j, k + 1) + x(i + 1, j, k + 1) + x(i, j + 1, k + 1) + x(i + 1, j + 1, k + 1)),
                            0.25 * (y(i, j, k + 1) + y(i + 1, j, k + 1) + y(i, j + 1, k + 1) + y(i + 1, j + 1, k + 1)),
                            0.25 * (z(i, j, k + 1) + z(i + 1, j, k + 1) + z(i, j + 1, k + 1) + z(i + 1, j + 1, k + 1))};
                        V += dot(Af_zp, minus(xf_zp, xc));

                        Jac(i, j, k, 0) = V / 3.0;
                    }
        }
    }

    // GCL Test
    // {
    //     auto &Jac_ = field("Jac");
    //     auto &Axi_ = field("JDxi");
    //     auto &Aeta_ = field("JDet");
    //     auto &Azeta_ = field("JDze");
    //     for (int ib = 0; ib < grd->nblock; ++ib)
    //     {
    //         auto &Jac = Jac_[ib];
    //         auto &A_xi = Axi_[ib];
    //         auto &A_et = Aeta_[ib];
    //         auto &A_ze = Azeta_[ib];
    //         Int3 lo = Jac.inner_lo();
    //         Int3 hi = Jac.inner_hi();
    //         for (int i = lo.i; i < hi.i; ++i)
    //             for (int j = lo.j; j < hi.j; ++j)
    //                 for (int k = lo.k; k < hi.k; ++k)
    //                 {
    //                     double error0 = A_xi(i + 1, j, k, 0) - A_xi(i, j, k, 0);
    //                     double error1 = A_xi(i + 1, j, k, 1) - A_xi(i, j, k, 1);
    //                     double error2 = A_xi(i + 1, j, k, 2) - A_xi(i, j, k, 2);
    //                     error0 += A_et(i, j + 1, k, 0) - A_et(i, j, k, 0);
    //                     error1 += A_et(i, j + 1, k, 1) - A_et(i, j, k, 1);
    //                     error2 += A_et(i, j + 1, k, 2) - A_et(i, j, k, 2);
    //                     error0 += A_ze(i, j, k + 1, 0) - A_ze(i, j, k, 0);
    //                     error1 += A_ze(i, j, k + 1, 1) - A_ze(i, j, k, 1);
    //                     error2 += A_ze(i, j, k + 1, 2) - A_ze(i, j, k, 2);
    //                     std::cout << sqrt(error0 * error0 + error1 * error1 + error2 * error2) << "\n";
    //                 }
    //     }
    // }
}