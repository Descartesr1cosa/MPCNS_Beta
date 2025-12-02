#include "Boundary.h"
#include "array"

void MHD_Boundary::apply_cell_patch_U(const TOPO::PhysicalPatch &patch, int32_t field_id)
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
    else
    {
        // 默认给一个简单的拷贝边界
        apply_cell_copy(U, patch);
    }
}

void MHD_Boundary::apply_cell_copy(FieldBlock &U, const TOPO::PhysicalPatch &patch)
{
    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.this_box_node.lo; // 含 ghost 的逻辑下界
    const Int3 hi = patch.this_box_node.hi; // 含 ghost 的逻辑上界（不含）
    const int ncomp = desc.ncomp;

    int cycle[3] = {patch.raw->cycle[0], patch.raw->cycle[1], patch.raw->cycle[2]};

    // 根据 direction 判断偏移量
    int shift = (patch.direction > 0) ? 1 : 0;

    for (int g = 1 - shift; g <= ngh - shift; ++g) // ghost 层
    {
        for (int i = lo.i + cycle[0] * g; i < hi.i + cycle[0] * g; i++)
            for (int j = lo.j + cycle[1] * g; j < hi.j + cycle[1] * g; j++)
                for (int k = lo.k + cycle[2] * g; k < hi.k + cycle[2] * g; k++)
                {
                    for (int m = 0; m < ncomp; ++m)
                        U(i, j, k, m) = U(i - cycle[0] * (g + shift), j - cycle[1] * (g + shift), k - cycle[2] * (g + shift), m);
                }
    }
}

void MHD_Boundary::apply_cell_farfield(FieldBlock &U, const TOPO::PhysicalPatch &patch)
{
    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.this_box_node.lo; // 含 ghost 的逻辑下界
    const Int3 hi = patch.this_box_node.hi; // 含 ghost 的逻辑上界（不含）
    const int ncomp = desc.ncomp;

    int cycle[3] = {patch.raw->cycle[0], patch.raw->cycle[1], patch.raw->cycle[2]};

    // 根据 direction 判断偏移量
    int shift = (patch.direction > 0) ? 1 : 0;

    for (int g = 1 - shift; g <= ngh - shift; ++g) // ghost 层
    {
        for (int i = lo.i + cycle[0] * g; i < hi.i + cycle[0] * g; i++)
            for (int j = lo.j + cycle[1] * g; j < hi.j + cycle[1] * g; j++)
                for (int k = lo.k + cycle[2] * g; k < hi.k + cycle[2] * g; k++)
                {
                    U(i, j, k, 0) = farfield_rho;
                    U(i, j, k, 1) = farfield_rho * farfield_u;
                    U(i, j, k, 2) = farfield_rho * farfield_v;
                    U(i, j, k, 3) = farfield_rho * farfield_w;
                    U(i, j, k, 4) = 0.5 * farfield_rho * (farfield_u * farfield_u + farfield_v * farfield_v + farfield_w * farfield_w) + farfield_p / (gamma_ - 1.0);
                    U(i, j, k, 4) += 0.5 * inver_MA2 * (farfield_bx * farfield_bx + farfield_by * farfield_by + farfield_bz * farfield_bz);
                }
    }
}

void MHD_Boundary::apply_cell_wall(FieldBlock &U, FieldBlock &Bcell, const TOPO::PhysicalPatch &patch)
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
    const Int3 lo = patch.this_box_node.lo; // 含 ghost 的逻辑下界
    const Int3 hi = patch.this_box_node.hi; // 含 ghost 的逻辑上界（不含）
    const int ncomp = desc.ncomp;

    int cycle[3] = {patch.raw->cycle[0], patch.raw->cycle[1], patch.raw->cycle[2]};

    // 根据 direction 判断偏移量
    int shift = (patch.direction > 0) ? 1 : 0;
    double sign = (patch.direction > 0) ? -1.0 : 1.0; // 使得face的面积矢量指向流体侧面

    // 面积矢量
    FieldBlock *face;
    if (cycle[0] != 0)
        face = &fld_->field("JDxi", patch.this_block);
    else if (cycle[1] != 0)
        face = &fld_->field("JDet", patch.this_block);
    else if (cycle[2] != 0)
        face = &fld_->field("JDze", patch.this_block);

    double p0, u0, v0, w0, rho0, Bx, By, Bz, Eb;
    // 壁面层
    for (int i = lo.i - cycle[0] * shift; i < hi.i - cycle[0] * shift; i++)
        for (int j = lo.j - cycle[1] * shift; j < hi.j - cycle[1] * shift; j++)
            for (int k = lo.k - cycle[2] * shift; k < hi.k - cycle[2] * shift; k++)
            {
                rho0 = U(i, j, k, 0);
                std::array<double, 3> vec_face = {(*face)(i + cycle[0] * shift, j + cycle[1] * shift, k + cycle[2] * shift, 0),
                                                  (*face)(i + cycle[0] * shift, j + cycle[1] * shift, k + cycle[2] * shift, 1),
                                                  (*face)(i + cycle[0] * shift, j + cycle[1] * shift, k + cycle[2] * shift, 2)};
                double norm = fmax(sqrt(dot(vec_face, vec_face)), 1e-15);
                vec_face = scalar(vec_face, sign / norm);
                // 此时vec_face永远指向流体侧边
                std::array<double, 3> vec_v = {U(i, j, k, 1), U(i, j, k, 2), U(i, j, k, 3)};
                vec_v = scalar(vec_v, 1.0 / rho0);

                u0 = vec_v[0];
                v0 = vec_v[1];
                w0 = vec_v[2];

                // 消去法向分量，保留切向
                // double dot_product = dot(vec_v, vec_face);
                // vec_v = minus(vec_v, scalar(vec_face, dot_product));
                vec_v = {0.0, 0.0, 0.0};

                Bx = Bcell(i, j, k, 0);
                By = Bcell(i, j, k, 1);
                Bz = Bcell(i, j, k, 2);
                Eb = 0.5 * inver_MA2 * (Bx * Bx + By * By + Bz * Bz);
                p0 = (U(i, j, k, 4) - 0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0) - Eb) * (gamma_ - 1.0);

                U(i, j, k, 1) = rho0 * vec_v[0];
                U(i, j, k, 2) = rho0 * vec_v[1];
                U(i, j, k, 3) = rho0 * vec_v[2];
                U(i, j, k, 4) = p0 / (gamma_ - 1.0) + 0.5 * rho0 * dot(vec_v, vec_v) + Eb;
            }

    for (int g = 1 - shift; g <= ngh - shift; ++g) // ghost 层
    {
        for (int i = lo.i + cycle[0] * g; i < hi.i + cycle[0] * g; i++)
            for (int j = lo.j + cycle[1] * g; j < hi.j + cycle[1] * g; j++)
                for (int k = lo.k + cycle[2] * g; k < hi.k + cycle[2] * g; k++)
                {
                    for (int m = 0; m < ncomp; ++m)
                        U(i, j, k, m) = U(i - cycle[0] * (g + shift), j - cycle[1] * (g + shift), k - cycle[2] * (g + shift), m);
                }
    }
}

void MHD_Boundary::apply_cell_pole(FieldBlock &U, const TOPO::PhysicalPatch &patch)
{
    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.this_box_node.lo; // 含 ghost 的逻辑下界
    const Int3 hi = patch.this_box_node.hi; // 含 ghost 的逻辑上界（不含）
    const int ncomp = desc.ncomp;

    int cycle[3] = {patch.raw->cycle[0], patch.raw->cycle[1], patch.raw->cycle[2]};

    // 根据 direction 判断偏移量
    int shift = (patch.direction > 0) ? -1 : 0;

    int mz = fld_->grd->grids(patch.this_block).mz;
    double count;
    std::vector<double> temp;
    temp.resize(ncomp);
    for (int i = lo.i + cycle[0] * shift; i < hi.i + cycle[0] * shift; i++)
        for (int j = lo.j + cycle[1] * shift; j < hi.j + cycle[1] * shift; j++)
        {
            for (int m = 0; m < ncomp; ++m)
                temp[m] = 0.0;
            count = 0.0;

            // Sum
            for (int k = 0; k < mz; k++)
            {
                for (int m = 0; m < ncomp; ++m)
                    temp[m] += U(i, j, k, m);
                count += 1.0;
            }

            // Average
            for (int m = 0; m < ncomp; ++m)
                temp[m] /= count;

            // Unified
            for (int k = 0; k < mz; k++)
                for (int m = 0; m < ncomp; ++m)
                    U(i, j, k, m) = temp[m];
        }

    apply_cell_copy(U, patch);
}

void MHD_Boundary::apply_face_patch_B(const TOPO::PhysicalPatch &patch, int32_t field_id, StaggerLocation location)
{
    const int ib = patch.this_block;
    FieldBlock &U = fld_->field(field_id, ib); // 该face上的 U

    int loc = -1;
    switch (location)
    {
    case StaggerLocation::FaceXi:
        loc = 0;
        break;
    case StaggerLocation::FaceEt:
        loc = 1;
        break;
    case StaggerLocation::FaceZe:
        loc = 2;
        break;
    default:
        break;
    }

    // 这里你可以根据 patch.bc_name 判断是什么类型的边界
    if (patch.bc_name == "Farfield")
    {
        apply_face_farfield(U, patch, loc);
    }
    else
    {
        // 默认给一个简单的拷贝边界
        apply_face_copy(U, patch, loc);
    }
}

void MHD_Boundary::apply_face_copy(FieldBlock &U, const TOPO::PhysicalPatch &patch, int loc)
{
    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.this_box_node.lo; // 含 ghost 的逻辑下界
    const Int3 hi = patch.this_box_node.hi; // 含 ghost 的逻辑上界（不含）
    const int ncomp = desc.ncomp;

    int cycle[3] = {patch.raw->cycle[0], patch.raw->cycle[1], patch.raw->cycle[2]};

    // 根据 direction 判断偏移量
    // int shift = (patch.direction > 0) ? 1 : 0;

    int shift = -10086;

    if (cycle[loc] != 0) // face方向与边界面法向一致，此方向范围[0,mx] 无需平移1 与node相同
        shift = 0;
    else // face方向与边界面法向不一致，此方向范围[0,mx) 大号面需平移1 与cell相同
        shift = (patch.direction > 0) ? 1 : 0;

    for (int g = 1 - shift; g <= ngh - shift; ++g) // ghost 层
    {
        for (int i = lo.i + cycle[0] * g; i < hi.i + cycle[0] * g; i++)
            for (int j = lo.j + cycle[1] * g; j < hi.j + cycle[1] * g; j++)
                for (int k = lo.k + cycle[2] * g; k < hi.k + cycle[2] * g; k++)
                {
                    for (int m = 0; m < ncomp; ++m)
                        U(i, j, k, m) = U(i - cycle[0] * (g + shift), j - cycle[1] * (g + shift), k - cycle[2] * (g + shift), m);
                }
    }
}

void MHD_Boundary::apply_face_farfield(FieldBlock &U, const TOPO::PhysicalPatch &patch, int loc)
{
    const FieldDescriptor &desc = U.descriptor();
    const int ngh = desc.nghost;
    const Int3 lo = patch.this_box_node.lo; // 含 ghost 的逻辑下界
    const Int3 hi = patch.this_box_node.hi; // 含 ghost 的逻辑上界（不含）
    const int ncomp = desc.ncomp;

    int cycle[3] = {patch.raw->cycle[0], patch.raw->cycle[1], patch.raw->cycle[2]};

    // 面积矢量
    FieldBlock *face;
    if (cycle[0] != 0)
        face = &fld_->field("JDxi", patch.this_block);
    else if (cycle[1] != 0)
        face = &fld_->field("JDet", patch.this_block);
    else if (cycle[2] != 0)
        face = &fld_->field("JDze", patch.this_block);

    int shift = -10086;

    if (cycle[loc] != 0) // face方向与边界面法向一致，此方向范围[0,mx] 无需平移1 与node相同
        shift = 0;
    else // face方向与边界面法向不一致，此方向范围[0,mx) 大号面需平移1 与cell相同
        shift = (patch.direction > 0) ? 1 : 0;

    for (int g = 1 - shift; g <= ngh - shift; ++g) // ghost 层
    {
        for (int i = lo.i + cycle[0] * g; i < hi.i + cycle[0] * g; i++)
            for (int j = lo.j + cycle[1] * g; j < hi.j + cycle[1] * g; j++)
                for (int k = lo.k + cycle[2] * g; k < hi.k + cycle[2] * g; k++)
                {
                    // 使用壁面的面积矢量
                    // std::array<double, 3> vec_face = {(*face)(i - cycle[0] * (g + shift), j - cycle[1] * (g + shift), k - cycle[2] * (g + shift), 0),
                    //   (*face)(i - cycle[0] * (g + shift), j - cycle[1] * (g + shift), k - cycle[2] * (g + shift), 1),
                    //   (*face)(i - cycle[0] * (g + shift), j - cycle[1] * (g + shift), k - cycle[2] * (g + shift), 2)};

                    // U(i, j, k, 0) = farfield_bx * vec_face[0] + farfield_by * vec_face[1] + farfield_bz * vec_face[2];
                    U(i, j, k, 0) = 0.0;
                }
    }
}
