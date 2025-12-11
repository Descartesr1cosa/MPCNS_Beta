#include "Boundary.h"
#include "array"

void MHD_Boundary::apply_cell_copy(FieldBlock &U, PhysicalRegion &patch)
{
    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.box_bound.lo;
    const Int3 hi = patch.box_bound.hi;
    const Int3 cyc = patch.cycle;
    const int ncomp = desc.ncomp;

    int32_t i, j, k;
    for (int32_t ii = lo.i; ii < hi.i; ii++)
        for (int32_t jj = lo.j; jj < hi.j; jj++)
            for (int32_t kk = lo.k; kk < hi.k; kk++)
            {
                // ii jj kk在边界面上循环，属于计算网格范围
                for (int32_t ng = 1; ng <= ngh; ng++)
                {
                    // i j k在虚网格上循环
                    i = ii + ng * cyc.i;
                    j = jj + ng * cyc.j;
                    k = kk + ng * cyc.k;

                    for (int m = 0; m < ncomp; ++m)
                        U(i, j, k, m) = U(ii, jj, kk, m);
                }
            }
}

void MHD_Boundary::apply_cell_farfield(FieldBlock &U, PhysicalRegion &patch)
{
    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.box_bound.lo;
    const Int3 hi = patch.box_bound.hi;
    const Int3 cyc = patch.cycle;
    const int ncomp = desc.ncomp;

    int32_t i, j, k;
    for (int32_t ii = lo.i; ii < hi.i; ii++)
        for (int32_t jj = lo.j; jj < hi.j; jj++)
            for (int32_t kk = lo.k; kk < hi.k; kk++)
            {
                // ii jj kk在边界面上循环，属于计算网格范围
                for (int32_t ng = 1; ng <= ngh; ng++)
                {
                    // i j k在虚网格上循环
                    i = ii + ng * cyc.i;
                    j = jj + ng * cyc.j;
                    k = kk + ng * cyc.k;

                    U(i, j, k, 0) = farfield_rho;
                    U(i, j, k, 1) = farfield_rho * farfield_u;
                    U(i, j, k, 2) = farfield_rho * farfield_v;
                    U(i, j, k, 3) = farfield_rho * farfield_w;
                    U(i, j, k, 4) = 0.5 * farfield_rho * (farfield_u * farfield_u + farfield_v * farfield_v + farfield_w * farfield_w) + farfield_p / (gamma_ - 1.0);
                    U(i, j, k, 4) += 0.5 * inver_MA2 * (farfield_bx * farfield_bx + farfield_by * farfield_by + farfield_bz * farfield_bz);
                }
            }
}

void MHD_Boundary::apply_cell_wall(FieldBlock &U, FieldBlock &Bcell, PhysicalRegion &patch)
{
    auto dot = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    auto scalar = [&](const std::array<double, 3> &a, double sca) -> std::array<double, 3>
    {
        std::array<double, 3> b = {a[0] * sca, a[1] * sca, a[2] * sca};
        return b;
    };
    auto minus = [&](const std::array<double, 3> &a, const std::array<double, 3> &b) -> std::array<double, 3>
    {
        std::array<double, 3> c = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
        return c;
    };

    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.box_bound.lo;
    const Int3 hi = patch.box_bound.hi;
    const Int3 cyc = patch.cycle;

    double sign = (patch.raw->direction > 0) ? -1.0 : 1.0; // 使得face的面积矢量指向流体侧面
    int shift = (patch.raw->direction > 0) ? 1 : 0;        // 根据 direction 判断偏移量，获得壁面法向矢量

    // 面积矢量
    FieldBlock *face = nullptr;
    if (cyc.i != 0)
        face = &fld_->field("JDxi", patch.this_block);
    else if (cyc.j != 0)
        face = &fld_->field("JDet", patch.this_block);
    else if (cyc.k != 0)
        face = &fld_->field("JDze", patch.this_block);

    // 小的压力下限，只防止浮点抖动出负压，不去“硬改物理”
    const double p_floor = 1e-8 * farfield_p;

    double p0, rho0, Bx, By, Bz, Eb, kin;
    double B_add_x = fld_->par->GetDou("B_add_x");
    double B_add_y = fld_->par->GetDou("B_add_y");
    double B_add_z = fld_->par->GetDou("B_add_z");

    int32_t i, j, k;
    for (int32_t ii = lo.i; ii < hi.i; ii++)
        for (int32_t jj = lo.j; jj < hi.j; jj++)
            for (int32_t kk = lo.k; kk < hi.k; kk++)
            {
                // 获得壁面法向矢量
                std::array<double, 3> vec_face = {(*face)(ii + cyc.i * shift, jj + cyc.j * shift, kk + cyc.k * shift, 0),
                                                  (*face)(ii + cyc.i * shift, jj + cyc.j * shift, kk + cyc.k * shift, 1),
                                                  (*face)(ii + cyc.i * shift, jj + cyc.j * shift, kk + cyc.k * shift, 2)};
                double norm = fmax(sqrt(dot(vec_face, vec_face)), 1e-15);
                // 调整vec_face方向，使得其永远指向流体侧边
                vec_face = scalar(vec_face, sign / norm);

                // 获取壁面法向计算区域第一层流体参数
                rho0 = U(ii, jj, kk, 0);
                std::array<double, 3> vec_v = {U(ii, jj, kk, 1) / rho0, U(ii, jj, kk, 2) / rho0, U(ii, jj, kk, 3) / rho0};

                Bx = Bcell(ii, jj, kk, 0);
                By = Bcell(ii, jj, kk, 1);
                Bz = Bcell(ii, jj, kk, 2);

                kin = 0.5 * rho0 * dot(vec_v, vec_v);
                Eb = 0.5 * inver_MA2 * (Bx * Bx + By * By + Bz * Bz);
                p0 = (U(ii, jj, kk, 4) - kin - Eb) * (gamma_ - 1.0);
                if (p0 < p_floor)
                    p0 = p_floor; // 仅做极弱保护

                rho0 = 0.01;
                Bx = B_add_x;
                By = B_add_y;
                Bz = B_add_z;
                Eb = 0.5 * inver_MA2 * (Bx * Bx + By * By + Bz * Bz);

                // 反射法向速度，切向保留
                double vn2 = 2.0 * dot(vec_v, vec_face); // 2 (v·n)
                std::array<double, 3> v_g = minus(vec_v, scalar(vec_face, vn2));

                // double vn = dot(vec_v, vec_face);
                // // vec_face 指向流体侧
                // std::array<double, 3> v_g = vec_v;
                // if (vn > 0.0)
                // {
                //     // 去掉“从月球向流体侧”的法向分量
                //     v_g = minus(vec_v, scalar(vec_face, vn));
                // }

                // ii jj kk在边界面上循环，属于计算网格范围
                for (int32_t ng = 1; ng <= ngh; ng++)
                {
                    // i j k在虚网格上循环
                    i = ii + ng * cyc.i;
                    j = jj + ng * cyc.j;
                    k = kk + ng * cyc.k;
                    //===================================================================
                    // const double damp = 0.8; // 0.8~1.0 之间试
                    U(i, j, k, 0) = rho0;
                    U(i, j, k, 1) = rho0 * v_g[0]; // damp * rho0 * v_g[0];
                    U(i, j, k, 2) = rho0 * v_g[1]; // damp * rho0 * v_g[1];
                    U(i, j, k, 3) = rho0 * v_g[2]; // damp * rho0 * v_g[2];
                    U(i, j, k, 4) = p0 / (gamma_ - 1.0) + 0.5 * rho0 * dot(v_g, v_g) + Eb;
                }
            }
}

void MHD_Boundary::apply_cell_pole(FieldBlock &U, PhysicalRegion &patch)
{
    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.box_bound.lo;
    const Int3 hi = patch.box_bound.hi;
    const Int3 cyc = patch.cycle;
    const int ncomp = desc.ncomp;

    double count;
    std::vector<double> temp;
    temp.resize(ncomp);
    int mz = fld_->grd->grids(patch.this_block).mz;

    for (int32_t ii = lo.i; ii < hi.i; ii++)
        for (int32_t jj = lo.j; jj < hi.j; jj++)
        {
            for (int m = 0; m < ncomp; ++m)
                temp[m] = 0.0;
            count = 0.0;

            // Sum
            for (int kk = 0; kk < mz; kk++)
            {
                for (int m = 0; m < ncomp; ++m)
                    temp[m] += U(ii, jj, kk, m);
                count += 1.0;
            }

            // Average
            for (int m = 0; m < ncomp; ++m)
                temp[m] /= count;

            // Unified
            for (int kk = 0; kk < mz; kk++)
                for (int m = 0; m < ncomp; ++m)
                    U(ii, jj, kk, m) = temp[m];
        }

    apply_cell_copy(U, patch);
}

void MHD_Boundary::apply_face_copy(FieldBlock &U, PhysicalRegion &patch)
{
    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.box_bound.lo;
    const Int3 hi = patch.box_bound.hi;
    const Int3 cyc = patch.cycle;
    const int ncomp = desc.ncomp;

    int32_t i, j, k;
    for (int32_t ii = lo.i; ii < hi.i; ii++)
        for (int32_t jj = lo.j; jj < hi.j; jj++)
            for (int32_t kk = lo.k; kk < hi.k; kk++)
            {
                // ii jj kk在边界面上循环，属于计算网格范围
                for (int32_t ng = 1; ng <= ngh; ng++)
                {
                    // i j k在虚网格上循环
                    i = ii + ng * cyc.i;
                    j = jj + ng * cyc.j;
                    k = kk + ng * cyc.k;

                    for (int m = 0; m < ncomp; ++m)
                        U(i, j, k, m) = U(ii, jj, kk, m);
                }
            }
}

void MHD_Boundary::apply_derived_cell_wall(FieldBlock &U, PhysicalRegion &patch)
{

    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.box_bound.lo;
    const Int3 hi = patch.box_bound.hi;
    const Int3 cyc = patch.cycle;
    const int ncomp = desc.ncomp;

    double B_add_x = fld_->par->GetDou("B_add_x");
    double B_add_y = fld_->par->GetDou("B_add_y");
    double B_add_z = fld_->par->GetDou("B_add_z");

    int32_t i, j, k;
    for (int32_t ii = lo.i; ii < hi.i; ii++)
        for (int32_t jj = lo.j; jj < hi.j; jj++)
            for (int32_t kk = lo.k; kk < hi.k; kk++)
            {
                // ii jj kk在边界面上循环，属于计算网格范围
                for (int32_t ng = 1; ng <= ngh; ng++)
                {
                    // i j k在虚网格上循环
                    i = ii + ng * cyc.i;
                    j = jj + ng * cyc.j;
                    k = kk + ng * cyc.k;

                    U(i, j, k, 0) = B_add_x;
                    U(i, j, k, 1) = B_add_y;
                    U(i, j, k, 2) = B_add_z;
                }
            }
}
