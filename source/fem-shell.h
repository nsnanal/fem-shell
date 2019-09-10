/*************************************************\
*     ___  ___          __       ___              *
*    /__  /__  /\/\  _ /_   /_/ /__  /   /        *
*   /    /__  /    \   __/ / / /__  /__ /__       *
*                                                 *
*  developed by Stephan Herb in his master-thesis *
*  June - December, 2015                          *
*  Stand-alone version                            *
\*************************************************/

#ifndef FEMSHELL_H
#define FEMSHELL_H

// C++ include files that we need
#include <algorithm>
#include <iostream>
#include <math.h>
#include <unordered_set>

// libMesh includes
//#include "libmesh/perf_log.h" // for performance logging
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/getpot.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/zero_function.h"

namespace ShellSolid
{
  // Bring in everything from the libMesh namespace
  using namespace libMesh;

  struct shellparam
  {
    shellparam(){};
    ~shellparam(){};
    bool read_parameters(int, char **);
    std::string in_filename;    // mesh file for import
    std::string force_filename; // force file
    std::string out_filename;   // output file name
    bool debug;
    Real nu;        // Poisson's ratio
    Real em;        // Young's modulus
    Real thickness; // Mesh thickness
    bool isOutfileSet;
  };

  class shellsolid
  {
  public:
    shellsolid(const SerialMesh &, const shellparam &);
    ~shellsolid(){};
    // function prototypes:
    void read_forcing();

    void make_constraints(std::map<boundary_id_type, unsigned int> &);

    void writeOutput();

    void set_solution(const std::vector<Number> &);

    inline const std::vector<Number> &get_solution() { return this->sols; };

    void run();

  private:
    void initMaterialMatrices();

    static void initElement(const Elem **elem,
                            DenseMatrix<Real> &transUV,
                            DenseMatrix<Real> &trafo,
                            DenseMatrix<Real> &dphi,
                            Real *area);

    static void calcPlane(EquationSystems &es,
                          ElemType type,
                          DenseMatrix<Real> &transUV,
                          DenseMatrix<Real> &dphi,
                          Real *area,
                          DenseMatrix<Real> &Ke_m);

    static void calcPlate(EquationSystems &es,
                          ElemType type,
                          DenseMatrix<Real> &dphi,
                          Real *area,
                          DenseMatrix<Real> &Ke_p);

    void
    B_plane_tri(Real *area, DenseMatrix<Real> &dphi, DenseMatrix<Real> &out);

    void B_plane_quad(Real *area,
                      DenseMatrix<Real> &dphi,
                      DenseMatrix<Real> &transUV,
                      Real qp_x,
                      Real qp_y,
                      DenseMatrix<Real> &out);

    static void evalBTri(EquationSystems &es,
                         DenseVector<Real> &C,
                         Real L1,
                         Real L2,
                         DenseMatrix<Real> &dphi,
                         DenseMatrix<Real> &out);

    static void evalBQuad(EquationSystems &es,
                          DenseMatrix<Real> &Hcoeffs,
                          Real xi,
                          Real eta,
                          DenseMatrix<Real> &Jinv,
                          DenseMatrix<Real> &out);

    static void constructStiffnessMatrix(EquationSystems &es,
                                         ElemType type,
                                         DenseMatrix<Real> &Ke_m,
                                         DenseMatrix<Real> &Ke_p,
                                         DenseMatrix<Real> &K_out);

    static void localToGlobalTrafo(EquationSystems &es,
                                   ElemType type,
                                   DenseMatrix<Real> &trafo,
                                   DenseMatrix<Real> &Ke_inout);

    static void contribRHS(EquationSystems &es,
                           const Elem **elem,
                           DenseVector<Real> &Fe,
                           std::unordered_set<unsigned int> *processedNodes);

    static void assemble_elasticity(EquationSystems &es,
                                    const std::string &system_name);

    void stress_calculation();

    SerialMesh mesh;
    EquationSystems equation_systems;
    LinearImplicitSystem &system;
    ExplicitSystem &stress_system;   // Membrane Stresses
    ExplicitSystem &stress_system_b; // Bending Stresses

    std::string in_filename;               // mesh file for import
    std::string force_filename;            // forcefile
    std::string out_filename;              // output file name
    bool debug;                            // show debug messages?
    Real nu;                               // Poisson's ratio
    Real em;                               // Young's modulus
    Real thickness;                        // Mesh thickness
    std::vector<DenseVector<Real>> forces; // nodal force vector
    std::vector<Number> sols;              // solution vector
    bool isOutfileSet;                     // should outputs be written or not

    DenseMatrix<Real> Dp, Dm; // material matrix for plate (Dp) and plane (Dm)
  };
} // namespace ShellSolid
#endif // FEMSHELL_H
