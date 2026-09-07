#ifndef AMREX_EB2_DISTRIBUTEDGEOMETRYSHOP_H_
#define AMREX_EB2_DISTRIBUTEDGEOMETRYSHOP_H_

#include <AMReX_Config.H>
#include <AMReX_EB2_IF_Base.H>
#include <AMReX_EB2_Graph.H>
#include <AMReX_Geometry.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Print.H>
#include <AMReX_Array.H>
#include <memory>
#include <type_traits>
#include <cmath>
#include <utility>
#include <vector>

namespace amrex::EB2 {

/**
 * \brief Functor that interpolates a Signed Distance Function from a multilevel 
 *        set of MultiFabs.
 */
class DistributedSDF {
public:
    DistributedSDF(Vector<MultiFab> mfs, Vector<Geometry> geoms)
        : m_mfs(std::move(mfs)), m_geoms(std::move(geoms)) {}

    AMREX_GPU_HOST_DEVICE Real operator()(Real x, Real y, Real z) const {
        int lev = -1;
        // Iterate from finest to coarsest to ensure we get the highest resolution data.
        for (int l = (int)m_mfs.size() - 1; l >= 0; --l) {
            // 1. Check if point is within the global probability domain of this level.
            if (m_geoms[l].ProbLo()[0] <= x && x <= m_geoms[l].ProbHi()[0] &&
                m_geoms[l].ProbLo()[1] <= y && y <= m_geoms[l].ProbHi()[1] &&
                m_geoms[l].ProbLo()[2] <= z && z <= m_geoms[l].ProbHi()[2]) {
                
                // 2. Verify that the point actually falls within one of the defined boxes (FArrayBox)
                // because the domain might have "holes" where data is not refined to this level.
                const auto& boxes = m_mfs[l].boxArray();
                bool found = false;
                
                // In a real sampling scenario, we map coordinates to a global index.
                Real dx = m_geoms[l].CellSize()[0];
                Real dy = m_geoms[l].CellSize()[1];
                Real dz = m_geoms[l].CellSize()[2];
                
                // Map physical coords to indices (cell-centered)
                int i = static_cast<int>(std::floor((x - m_geoms[l].ProbLo()[0]) / dx - 0.5_rt));
                int j = static_cast<int>(std::floor((y - m_geoms[l].ProbLo()[1]) / dy - 0.5_rt));
                int k = static_cast<int>(std::floor((z - m_geoms[l].ProbLo()[2]) / dz - 0.5_rt));
                IntVect iv(i, j, k);

                for (int b = 0; b < (int)boxes.size(); ++b) {
                    if (boxes[b].contains(iv)) {
                        found = true;
                        break;
                    }
                }
                
                if (found) {
                    lev = l;
                    break;
                }
            }
        }

        if (lev == -1) return 1.0_rt; // Outside all defined SDF regions

        const auto& geom = m_geoms[lev];
        const auto& mf = m_mfs[lev];
        
        Real dx = geom.CellSize()[0];
        Real dy = geom.CellSize()[1];
        Real dz = geom.CellSize()[2];
        
        Real fx = (x - geom.ProbLo()[0]) / dx - 0.5_rt;
        Real fy = (y - geom.ProbLo()[1]) / dy - 0.5_rt;
        Real fz = (z - geom.ProbLo()[2]) / dz - 0.5_rt;

        int i = static_cast<int>(std::floor(fx));
        int j = static_cast<int>(std::floor(fy));
        int k = static_cast<int>(std::floor(fz));

        Real tx = fx - i;
        Real ty = fy - j;
        Real tz = fz - k;

        // Helper to find the correct FArrayBox for a given global index.
        auto get_val = [&](int ii, int jj, int kk) {
            IntVect iv(ii, jj, kk);
            const auto& boxes = mf.boxArray();
            for (int b = 0; b < (int)boxes.size(); ++b) {
                if (boxes[b].contains(iv)) {
                    return mf[b](iv);
                }
            }
            // If out of bounds of current level's boxes, we fall back to 0 or 
            // technically we should have stayed at a coarser level.
            return 0.0_rt; 
        };

        Real v000 = get_val(i, j, k);
        Real v100 = get_val(i+1, j, k);
        Real v010 = get_val(i, j+1, k);
        Real v110 = get_val(i+1, j+1, k);
        Real v001 = get_val(i, j, k+1);
        Real v101 = get_val(i+1, j, k+1);
        Real v011 = get_val(i, j+1, k+1);
        Real v111 = get_val(i+1, j+1, k+1);

        return (1-tx)*(1-ty)*(1-tz)*v000 + tx*(1-ty)*(1-tz)*v100 +
               (1-tx)*ty*(1-tz)*v010 + tx*ty*(1-tz)*v110 +
               (1-tx)*(1-ty)*tz*v001 + tx*(1-ty)*tz*v101 +
               (1-tx)*ty*tz*v011 + tx*ty*tz*v111;
    }

private:
    Vector<MultiFab> m_mfs;
    Vector<Geometry> m_geoms;
};

template <class F>
AMREX_GPU_HOST_DEVICE
Real
IF_f_dist (F const& f, GpuArray<Real,AMREX_SPACEDIM> const& p)
{
    return f(AMREX_D_DECL(p[0],p[1],p[2]));
}

template <class F>
AMREX_GPU_HOST_DEVICE
Real
BrentRootFinder_dist (GpuArray<Real,AMREX_SPACEDIM> const& lo,
                       GpuArray<Real,AMREX_SPACEDIM> const& hi,
                       int rangedir, F const& f)
{
    const Real tol = 1.e-12;
    const Real EPS = 3.0e-15;
    const int MAXITER = 100;

    GpuArray<Real,AMREX_SPACEDIM> aPt = lo;
    GpuArray<Real,AMREX_SPACEDIM> bPt = hi;

    Real fa = IF_f_dist(f, aPt);
    Real fb = IF_f_dist(f, bPt);
    Real c = bPt[rangedir];
    Real fc = fb;

    if (fb*fa > 0.0_rt) return 0.0_rt;
    if (fa == 0.0_rt) return aPt[rangedir];
    if (fb == 0.0_rt) return bPt[rangedir];

    Real d = 0.0_rt, e = 0.0_rt;
    for (int i = 0; i < MAXITER; ++i) {
        if (fb*fc > 0) {
            c = aPt[rangedir]; fc = fa;
            d = bPt[rangedir] - aPt[rangedir]; e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
            aPt[rangedir] = bPt[rangedir]; bPt[rangedir] = c; c = aPt[rangedir];
            fa = fb; fb = fc; fc = fa;
        }
        Real tol1 = 2.0_rt * EPS * std::abs(bPt[rangedir]) + 0.5_rt * tol;
        Real xm = 0.5_rt * (c - bPt[rangedir]);
        if (std::abs(xm) <= tol1 || fb == 0.0_rt) break;
        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
            Real s = fb / fa;
            Real p, q, r;
            if (aPt[rangedir] == c) { p = 2.0_rt * xm * s; q = 1.0_rt - s; }
            else {
                q = fa / fc; r = fb / fc;
                p = s * (2.0_rt * xm * q * (q-r) - (bPt[rangedir]-aPt[rangedir]) * (r-1.0_rt));
                q = (q-1.0_rt) * (r-1.0_rt) * (s-1.0_rt);
            }
            if (p > 0) q = -q;
            p = std::abs(p);
            if (2.0_rt * p < amrex::min(3.0_rt*xm*q-std::abs(tol1*q), 1.0_rt*std::abs(e*q))) {
                e = d; d = p / q;
            } else { d = xm; e = d; }
        } else { d = xm; e = d; }
        aPt[rangedir] = bPt[rangedir]; fa = fb;
        if (std::abs(d) > tol1) bPt[rangedir] = bPt[rangedir] + d;
        else bPt[rangedir] = (xm < 0) ? bPt[rangedir] - tol1 : bPt[rangedir] + tol1;
        fb = IF_f_dist(f, bPt);
    }
    return bPt[rangedir];
}

template <class F = DistributedSDF>
class DistributedGeometryShop {
public:
    static constexpr int in_fluid = -1;
    static constexpr int on_boundary = 0;
    static constexpr int in_body = 1;
    static constexpr int allregular = -1;
    static constexpr int mixedcells = 0;
    static constexpr int allcovered = 1;

    using FunctionType = F;

    explicit DistributedGeometryShop(F f) : m_f(std::move(f)) {}

    [[nodiscard]] F const& GetImpFunc() const& { return m_f; }

    [[nodiscard]] int getBoxType_Cpu(const Box& bx, Geometry const& geom) const noexcept {
        return DistributedGeometryShop::mixedcells;
    }

    void getIntercept_Cpu(Array4<Real> const& inter,
                          Array4<Type_t const> const& type,
                          Geometry const& geom,
                          Box const& bounding_box,
                          const int idim) const noexcept {
        auto const& dx = geom.CellSizeArray();
        auto const& problo = geom.ProbLoArray();
        const auto& blo = amrex::lbound(bounding_box);
        const auto& bhi = amrex::ubound(bounding_box);
        const Box bx{inter};
        
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept {
            if (type(i,j,k) == Type::irregular) {
                IntVect ivlo(AMREX_D_DECL(i,j,k));
                IntVect ivhi(AMREX_D_DECL(i,j,k));
                ivhi[idim] += 1;
                inter(i,j,k) = BrentRootFinder_dist(
                    {AMREX_D_DECL(problo[0]+amrex::Clamp(ivlo[0],blo.x,bhi.x)*dx[0],
                                   problo[1]+amrex::Clamp(ivlo[1],blo.y,bhi.y)*dx[1],
                                   problo[2]+amrex::Clamp(ivlo[2],blo.z,bhi.z)*dx[2])},
                    {AMREX_D_DECL(problo[0]+amrex::Clamp(ivhi[0],blo.x,bhi.x)*dx[0],
                                   problo[1]+amrex::Clamp(ivhi[1],blo.y,bhi.y)*dx[1],
                                   problo[2]+amrex::Clamp(ivhi[2],blo.z,bhi.z)*dx[2])},
                    idim, m_f);
            } else {
                inter(i,j,k) = std::numeric_limits<Real>::quiet_NaN();
            }
        });
    }

private:
    F m_f;
};

template <class F = DistributedSDF>
DistributedGeometryShop<F> makeDistributedShop(F&& f) {
    return DistributedGeometryShop<std::decay_t<F>>(std::forward<F>(f));
}

} // namespace amrex::EB2

#endif
