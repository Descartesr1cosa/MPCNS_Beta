#include "Boundary.h"
#include "array"

void MHD_Boundary::apply_cell_patch_U(PhysicalRegion &patch, int32_t field_id)
{
    const int ib = patch.this_block;
    FieldBlock &U = fld_->field(field_id, ib);      // 该块上的 U
    FieldBlock &B_cell = fld_->field("B_cell", ib); // 该块上的 U

    // 这里你可以根据 patch.bc_name 判断是什么类型的边界
    if (patch.bc_name == "Solid_Surface")
    {
        apply_cell_wall(U, B_cell, patch);
    }
    else if (patch.bc_name == "Outflow")
    {
        // apply_bc_outflow(U, patch);
        apply_cell_copy(U, patch);
    }
    else if (patch.bc_name == "Farfield")
    {
        apply_cell_farfield(U, patch);
    }
    else if (patch.bc_name == "Pole")
        apply_cell_pole(U, patch);
    else
    {
        // 默认给一个简单的拷贝边界
        apply_cell_copy(U, patch);
    }
}

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
                for (int32_t ng = 1; ng <= ngg; ng++)
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
                for (int32_t ng = 1; ng <= ngg; ng++)
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
    const int ncomp = desc.ncomp;

    double sign = (patch.raw->direction > 0) ? -1.0 : 1.0; // 使得face的面积矢量指向流体侧面

    int shift = (patch.raw->direction > 0) ? 1 : 0; // 根据 direction 判断偏移量，获得壁面法向矢量
    // 面积矢量
    FieldBlock *face;
    if (cyc.i != 0)
        face = &fld_->field("JDxi", patch.this_block);
    else if (cyc.j != 0)
        face = &fld_->field("JDet", patch.this_block);
    else if (cyc.k != 0)
        face = &fld_->field("JDze", patch.this_block);

    double p0, u0, v0, w0, rho0, Bx, By, Bz, Eb;
    double B_add_x = fld_->par->GetDou("B_add_x");
    double B_add_y = fld_->par->GetDou("B_add_y");
    double B_add_z = fld_->par->GetDou("B_add_z");

    int32_t i, j, k;
    for (int32_t ii = lo.i; ii < hi.i; ii++)
        for (int32_t jj = lo.j; jj < hi.j; jj++)
            for (int32_t kk = lo.k; kk < hi.k; kk++)
            {
                // ii jj kk在边界面上循环，属于计算网格范围
                for (int32_t ng = 1; ng <= ngg; ng++)
                {
                    // i j k在虚网格上循环
                    i = ii + ng * cyc.i;
                    j = jj + ng * cyc.j;
                    k = kk + ng * cyc.k;

                    //===================================================================
                    // 获得壁面法向矢量
                    std::array<double, 3> vec_face = {(*face)(ii + cyc.i * shift, jj + cyc.j * shift, kk + cyc.k * shift, 0),
                                                      (*face)(ii + cyc.i * shift, jj + cyc.j * shift, kk + cyc.k * shift, 1),
                                                      (*face)(ii + cyc.i * shift, jj + cyc.j * shift, kk + cyc.k * shift, 2)};
                    double norm = fmax(sqrt(dot(vec_face, vec_face)), 1e-15);
                    // 调整vec_face方向，使得其永远指向流体侧边
                    vec_face = scalar(vec_face, sign / norm);

                    // 获取壁面法向计算区域第一层流体参数
                    rho0 = U(ii, jj, kk, 0);
                    std::array<double, 3> vec_v = {U(ii, jj, kk, 1), U(ii, jj, kk, 2), U(ii, jj, kk, 3)};
                    vec_v = scalar(vec_v, 1.0 / rho0);
                    u0 = vec_v[0];
                    v0 = vec_v[1];
                    w0 = vec_v[2];
                    Bx = Bcell(ii, jj, kk, 0);
                    By = Bcell(ii, jj, kk, 1);
                    Bz = Bcell(ii, jj, kk, 2);
                    Eb = 0.5 * inver_MA2 * (Bx * Bx + By * By + Bz * Bz);
                    p0 = (U(ii, jj, kk, 4) - 0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0) - Eb) * (gamma_ - 1.0);

                    p0 = fmax(p0, 0.001); // floor防止压力过低，使用无压力梯度
                    rho0 = 0.1;           // 吸收边界
                    // 虚网格法向分量反射，保留切向
                    double dot_product = 2.0 * dot(vec_v, vec_face);
                    vec_v = minus(vec_v, scalar(vec_face, dot_product));
                    // 壁面虚网格磁场直接忽略感应部分
                    Bx = B_add_x;
                    By = B_add_y;
                    Bz = B_add_z;
                    Eb = 0.5 * inver_MA2 * (Bx * Bx + By * By + Bz * Bz);

                    U(i, j, k, 0) = rho0;
                    U(i, j, k, 1) = rho0 * vec_v[0];
                    U(i, j, k, 2) = rho0 * vec_v[1];
                    U(i, j, k, 3) = rho0 * vec_v[2];
                    U(i, j, k, 4) = p0 / (gamma_ - 1.0) + 0.5 * rho0 * dot(vec_v, vec_v) + Eb;
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

void MHD_Boundary::apply_face_patch_B(PhysicalRegion &patch, int32_t field_id, StaggerLocation location)
{
    const int ib = patch.this_block;
    FieldBlock &U = fld_->field(field_id, ib); // 该face上的 U

    // int loc = -1;
    // switch (location)
    // {
    // case StaggerLocation::FaceXi:
    //     loc = 0;
    //     break;
    // case StaggerLocation::FaceEt:
    //     loc = 1;
    //     break;
    // case StaggerLocation::FaceZe:
    //     loc = 2;
    //     break;
    // default:
    //     break;
    // }

    // 默认给一个简单的拷贝边界
    apply_face_copy(U, patch);
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
                for (int32_t ng = 1; ng <= ngg; ng++)
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

void MHD_Boundary::apply_cell_wall_B(FieldBlock &U, PhysicalRegion &patch)
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
                for (int32_t ng = 1; ng <= ngg; ng++)
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
