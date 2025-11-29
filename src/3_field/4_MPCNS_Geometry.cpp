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
                        Area = minus(cross(rx, ry), cross(rx_, ry_));
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
                        Area = minus(cross(rx, ry), cross(rx_, ry_));
                        data(i, j, k, 0) = 0.5 * Area[0];
                        data(i, j, k, 1) = 0.5 * Area[1];
                        data(i, j, k, 2) = 0.5 * Area[2];
                    }
        }
    }

    {
        // Calc_ JDze  = Jac\nabla\eta = Area
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
                        Area = minus(cross(rx, ry), cross(rx_, ry_));
                        data(i, j, k, 0) = 0.5 * Area[0];
                        data(i, j, k, 1) = 0.5 * Area[1];
                        data(i, j, k, 2) = 0.5 * Area[2];
                    }
        }
    }

    {
        // Calc_ Jac = \SUM Area \cdot dr
        auto &data_ = field("Jac");
        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &data = data_[ib];

            auto &A_xi = field("JDxi")[ib];
            auto &A_eta = field("JDet")[ib];
            auto &A_zeta = field("JDze")[ib];

            auto &dual_x = grd->grids(ib).dual_x;
            auto &dual_y = grd->grids(ib).dual_y;
            auto &dual_z = grd->grids(ib).dual_z;

            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            Int3 inner_range_lo = data.inner_lo();
            Int3 inner_range_hi = data.inner_hi();
            std::array<double, 3> Area, center, face_point, dr;
            double Jacobian = 0.0;
            for (int i = inner_range_lo.i; i < inner_range_hi.i; i++)
                for (int j = inner_range_lo.j; j < inner_range_hi.j; j++)
                    for (int k = inner_range_lo.k; k < inner_range_hi.k; k++)
                    {
                        center = {dual_x(i + 1, j + 1, k + 1), dual_y(i + 1, j + 1, k + 1), dual_z(i + 1, j + 1, k + 1)};
                        Jacobian = 0.0;

                        // Minus
                        face_point = {x(i, j, k), y(i, j, k), z(i, j, k)};
                        dr = minus(center, face_point);

                        // Face_Xi minus
                        Area = {A_xi(i, j, k, 0), A_xi(i, j, k, 1), A_xi(i, j, k, 2)};
                        Jacobian += dot(Area, dr);

                        // Face_Eta minus
                        Area = {A_eta(i, j, k, 0), A_eta(i, j, k, 1), A_eta(i, j, k, 2)};
                        Jacobian += dot(Area, dr);

                        // Face_Xi minus
                        Area = {A_zeta(i, j, k, 0), A_zeta(i, j, k, 1), A_zeta(i, j, k, 2)};
                        Jacobian += dot(Area, dr);

                        // Plus
                        face_point = {x(i + 1, j + 1, k + 1), x(i + 1, j + 1, k + 1), x(i + 1, j + 1, k + 1)};
                        dr = minus(center, face_point);

                        // Face_Xi plus
                        Area = {A_xi(i + 1, j, k, 0), A_xi(i + 1, j, k, 1), A_xi(i + 1, j, k, 2)};
                        Jacobian += dot(Area, dr);

                        // Face_Eta plus
                        Area = {A_eta(i, j + 1, k, 0), A_eta(i, j + 1, k, 1), A_eta(i, j + 1, k, 2)};
                        Jacobian += dot(Area, dr);

                        // Face_Zeta plus
                        Area = {A_zeta(i, j, k + 1, 0), A_zeta(i, j, k + 1, 1), A_zeta(i, j, k + 1, 2)};
                        Jacobian += dot(Area, dr);

                        data(i, j, k, 0) = Jacobian / 3.0;
                    }
        }
    }
}