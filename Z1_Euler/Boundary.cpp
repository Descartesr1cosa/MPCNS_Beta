#include "Boundary.h"
#include "array"

void Euler_Boundary::apply_cell_patch_U(const TOPO::PhysicalPatch &patch, int32_t field_id)
{
    const int ib = patch.this_block;
    FieldBlock &U = fld_->field(field_id, ib); // 该块上的 U

    // 这里你可以根据 patch.bc_name 判断是什么类型的边界
    if (patch.bc_name == "Solid_Surface")
    {
        apply_cell_wall(U, patch);
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

void Euler_Boundary::apply_cell_copy(FieldBlock &U, const TOPO::PhysicalPatch &patch)
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

void Euler_Boundary::apply_cell_farfield(FieldBlock &U, const TOPO::PhysicalPatch &patch)
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
                }
    }
}

void Euler_Boundary::apply_cell_wall(FieldBlock &U, const TOPO::PhysicalPatch &patch)
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

    double p_inner, rho_inner, u_inner, v_inner, w_inner, rho0;
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

                // 消去法向分量，保留切向
                double dot_product = dot(vec_v, vec_face);
                vec_v = minus(vec_v, scalar(vec_face, dot_product));

                // rho_inner = U(i - cycle[0], j - cycle[1], k - cycle[2], 0);
                // u_inner = U(i - cycle[0], j - cycle[1], k - cycle[2], 1) / rho_inner;
                // v_inner = U(i - cycle[0], j - cycle[1], k - cycle[2], 2) / rho_inner;
                // w_inner = U(i - cycle[0], j - cycle[1], k - cycle[2], 3) / rho_inner;
                // p_inner = (U(i - cycle[0], j - cycle[1], k - cycle[2], 4) - 0.5 * rho_inner * (u_inner * u_inner + v_inner * v_inner + w_inner * w_inner)) * (gamma_ - 1.0);

                rho_inner = U(i, j, k, 0);
                u_inner = U(i, j, k, 1) / rho_inner;
                v_inner = U(i, j, k, 2) / rho_inner;
                w_inner = U(i, j, k, 3) / rho_inner;
                p_inner = (U(i, j, k, 4) - 0.5 * rho_inner * (u_inner * u_inner + v_inner * v_inner + w_inner * w_inner)) * (gamma_ - 1.0);

                U(i, j, k, 1) = rho0 * vec_v[0];
                U(i, j, k, 2) = rho0 * vec_v[1];
                U(i, j, k, 3) = rho0 * vec_v[2];
                U(i, j, k, 4) = p_inner / (gamma_ - 1.0) + 0.5 * rho0 * dot(vec_v, vec_v);
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
