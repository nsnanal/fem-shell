/*************************************************\
*     ___  ___          __       ___              *
*    /__  /__  /\/\  _ /_   /_/ /__  /   /        *
*   /    /__  /    \   __/ / / /__  /__ /__       *
*                                                 *
*  developed by Stephan Herb in his master-thesis *
*  June - December, 2015                          *
*  Stand-alone version                            *
\*************************************************/

#include "fem-shell.h"

namespace ShellSolid
{
  using namespace libMesh;
  /**
   * Processes the command-line arguments and fills the corresponding global
   * variables with values. If invalid arguments were set or too few, the
   * program gives out an error and quits.
   * @param argc Number of command-line arguments
   * @param argv Array storing the contents of the arguments
   * @return returns true if all parameters could be processes, false otherwise
   */
  bool shellparam::read_parameters(int argc, char **argv)
  {
    if (argc < 5)
      {
        err << "Error, must choose valid parameters.\n"
            << "Usage: " << argv[0] << " -nu -e -t -mesh [-out] [-d]\n"
            << "-nu:\t Possion's ratio (required)\n"
            << "-e:\t Elastic/Young's modulus E (required)\n"
            << "-t:\t Thickness (required)\n"
            << "-mesh:\t Input mesh file (*.xda/*.xdr or *.msh, required)\n"
            << "-out:\t Output file name (without extension, optional)\n"
            << "-d:\t Additional (debug) messages (1=on, 0=off (default))\n";

        return false;
      }

    // Parse command line
    bool failed = false;
    GetPot command_line(argc, argv);

    if (command_line.search(1, "-d"))
      debug = (command_line.next(0) == 1 ? true : false);

    if (command_line.search(1, "-nu"))
      nu = command_line.next(0.3);
    else
      {
        err << "ERROR: Poisson's ratio nu not specified!\n";
        failed = true;
      }

    if (command_line.search(1, "-e"))
      em = command_line.next(1.0e6);
    else
      {
        err << "ERROR: Elastic modulus E not specified!\n";
        failed = true;
      }

    if (command_line.search(1, "-t"))
      thickness = command_line.next(1.0);
    else
      {
        err << "ERROR: Mesh thickness t not specified!\n";
        failed = true;
      }

    if (command_line.search(1, "-mesh"))
      in_filename = command_line.next("mesh.xda");
    else
      {
        err << "ERROR: Mesh file not specified!\n";
        failed = true;
      }

    if (command_line.search(1, "-out"))
      {
        out_filename = command_line.next("out");
        isOutfileSet = true;
      }
    else
      isOutfileSet = false;

    // CONVENTION: force file has same file name as mesh file
    //             without extension, but with "_f" added at the end
    force_filename = in_filename;
    if (force_filename.find(".xda") != std::string::npos ||
        force_filename.find(".xdr") != std::string::npos ||
        force_filename.find(".msh") != std::string::npos)
      force_filename.resize(force_filename.size() - 4);

    force_filename += "_f";

    std::cout << "Run program with parameters:"
              << " debug messages = " << (debug ? "true" : "false")
              << ", nu = " << nu << ", E = " << em << ", t = " << thickness
              << ", mesh file = " << in_filename;
    if (isOutfileSet)
      std::cout << ", out-file = " << out_filename;
    std::cout << std::endl;
    return (!failed);
  }

  shellsolid::shellsolid(const SerialMesh &mesh, const shellparam &param)
    : mesh(mesh),
      equation_systems(this->mesh),
      system(equation_systems.add_system<LinearImplicitSystem>("Elasticity")),
      stress_system(
        equation_systems.add_system<ExplicitSystem>("StressSystem")),
      stress_system_b(
        equation_systems.add_system<ExplicitSystem>("StressSystem_b")),
      in_filename(param.in_filename),
      force_filename(param.force_filename),
      out_filename(param.out_filename),
      debug(param.debug),
      nu(param.nu),
      em(param.em),
      thickness(param.thickness),
      isOutfileSet(param.isOutfileSet)
  {
    equation_systems.parameters.insert<Real>("Thickness");
    equation_systems.parameters.set<Real>("Thickness") = thickness;
    equation_systems.parameters.insert<bool>("Debug");
    equation_systems.parameters.set<bool>("Debug") = debug;
  }

  void shellsolid::read_forcing()
  { // Load force file containing single nodal forces and moments
    std::filebuf fb;
    if (fb.open(force_filename.c_str(), std::ios::in))
      {
        std::istream input(&fb);
        int n_Forces;
        input >>
          n_Forces; // size of force vector (must equal number of mesh nodes)
        double factor = 1.0;
        input >>
          factor; // global factor multiplied to every single force and moment
        for (int i = 0; i < n_Forces; i++)
          {
            DenseVector<Real> p(6); // F_x, F_y, F_z, M_x, M_y, M_z
            for (int j = 0; j < 6; j++)
              input >> p(j);
            p *= factor;
            forces.push_back(p);
          }
      }
    equation_systems.parameters.insert<std::vector<DenseVector<Real>> *>(
      "Forcing terms");
    equation_systems.parameters.set<std::vector<DenseVector<Real>> *>(
      "Forcing terms") = &forces;
  }

  void
  shellsolid::make_constraints(std::map<boundary_id_type, unsigned int> &bcs)
  { // Add three displacement variables, u, v and w,
    // as well as three drilling variables theta_x, theta_y and theta_z to the
    // system
    unsigned int u_var = system.add_variable("u", FIRST, LAGRANGE);
    unsigned int v_var = system.add_variable("v", FIRST, LAGRANGE);
    unsigned int w_var = system.add_variable("w", FIRST, LAGRANGE);
    unsigned int tx_var = system.add_variable("tx", FIRST, LAGRANGE);
    unsigned int ty_var = system.add_variable("ty", FIRST, LAGRANGE);
    unsigned int tz_var = system.add_variable("tz", FIRST, LAGRANGE);

    // Add membrane stress tensor to the system
    stress_system.add_variable("sigma_xx", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_xy", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_xz", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_yx", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_yy", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_yz", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_zx", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_zy", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_zz", CONSTANT, MONOMIAL);

    // Add bending stress tensor to the system
    stress_system_b.add_variable("sigma_b_xx", CONSTANT, MONOMIAL);
    stress_system_b.add_variable("sigma_b_xy", CONSTANT, MONOMIAL);
    stress_system_b.add_variable("sigma_b_xz", CONSTANT, MONOMIAL);
    stress_system_b.add_variable("sigma_b_yx", CONSTANT, MONOMIAL);
    stress_system_b.add_variable("sigma_b_yy", CONSTANT, MONOMIAL);
    stress_system_b.add_variable("sigma_b_yz", CONSTANT, MONOMIAL);
    stress_system_b.add_variable("sigma_b_zx", CONSTANT, MONOMIAL);
    stress_system_b.add_variable("sigma_b_zy", CONSTANT, MONOMIAL);
    stress_system_b.add_variable("sigma_b_zz", CONSTANT, MONOMIAL);

    // We must add the Dirichlet boundary condition _before_ we call
    // equation_systems.init()
    // Create a ZeroFunction to initialize dirichlet_bc
    ConstFunction<Number> cf(0.0);
    for (auto &bc : bcs)
      {
        auto id = bc.first;
        auto flag = bc.second;
        std::vector<unsigned int> components;
        // 1-x, 2-y, 3-xy, 4-z, 5-xz, 6-yz, 7-xyz
        if (flag == 1 || flag == 3 || flag == 5 || flag == 7)
          {
            components.push_back(u_var);
          }
        if (flag == 2 || flag == 3 || flag == 6 || flag == 7)
          {
            components.push_back(v_var);
          }
        if (flag == 4 || flag == 5 || flag == 6 || flag == 7)
          {
            components.push_back(w_var);
          }
        std::set<boundary_id_type> tmp_id({id});
        DirichletBoundary dirichlet_bc(tmp_id, components, &cf);
        system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
      }
  }

  void shellsolid::set_solution(const std::vector<Number> &solutions)
  {
    const auto &dof_map = system.get_dof_map();
    std::vector<dof_id_type> dof_id_list(system.n_vars());
    for (auto node = mesh.active_nodes_begin(); node != mesh.active_nodes_end();
         ++node)
      {
        auto node_id = (*node)->id();
        dof_map.dof_indices(*node, dof_id_list);
        for (auto i : {0, 1, 2})
          {
            system.solution->set(dof_id_list[i], solutions[6 * node_id + i]);
          }
      }
  }

  void shellsolid::run()
  {
    initMaterialMatrices();

    read_forcing();

    // Initialize the data structures for the equation system.
    equation_systems.init();

    /**
     * Solve the system
     **/
    system.attach_assemble_function(shellsolid::assemble_elasticity);
    equation_systems.get_system("Elasticity").solve();

    equation_systems.build_solution_vector(sols);

    stress_calculation();

    if (debug)
      {
        std::cout << "System matrix:" << std::endl;
        system.matrix->print(std::cout);
        std::cout << std::endl << "RHS:" << std::endl;
        system.rhs->print(std::cout);
        std::cout << std::endl;
        // be sure that only the master process (id = 0) works on the solution,
        // since the rest of the processes only see their own partial solution
        if (global_processor_id() == 0)
          {
            std::cout << "Solution: u_vec = [";
            MeshBase::const_node_iterator no = mesh.nodes_begin();
            const MeshBase::const_node_iterator end_no = mesh.nodes_end();
            for (; no != end_no; ++no)
              {
                Node *nd = *no;
                int id = nd->id();
                Real displ_x = sols[24 * id];
                Real displ_y = sols[24 * id + 1];
                Real displ_z = sols[24 * id + 2];
                std::cout << "u= " << displ_x << ", v= " << displ_y
                          << ", w= " << displ_z;
                Real twist_x = sols[24 * id + 3];
                Real twist_y = sols[24 * id + 4];
                Real twist_z = sols[24 * id + 5];
                std::cout << ", tx= " << twist_x << ", ty= " << twist_y
                          << ", tz= " << twist_z << "]" << std::endl;
                // apply displacements to mesh nodes
                (*nd)(0) += displ_x;
                (*nd)(1) += displ_y;
                (*nd)(2) += displ_z;
              }
            std::cout << "]" << std::endl << std::endl;
          }
      }

    if (isOutfileSet)
      writeOutput();
  }

  /**
   * Initializes the material matrix for the plane element(s) (Dm) and plate
   * element(s) (Dp). Dp and Dm are global variables; nu,em,thickness are filled
   * with command-line argument values.
   */
  void shellsolid::initMaterialMatrices()
  {
    /*     /                   \
     *     | 1    nu      0    |
     * D = | nu   1       0    |
     *     | 0    0   (1-nu)/2 |
     *     \                   /
     */
    Dp.resize(3, 3);
    Dp(0, 0) = 1.0;
    Dp(0, 1) = nu;
    Dp(1, 0) = nu;
    Dp(1, 1) = 1.0;
    Dp(2, 2) = (1.0 - nu) / 2.0;
    Dm = Dp; // base matrix is same for Dm and Dp
    //         E
    // Dm = ------- * D
    //       1-nu²
    Dm *= em / (1.0 - nu * nu); // material matrix for plane part
    //          E * t³
    // Dp = ------------- * D
    //       12 * (1-nu²)
    Dp *= em * pow(thickness, 3.0) /
          (12.0 * (1.0 - nu * nu)); // material matrix for plate part
    // Insert the matrices into the equation system.
    equation_systems.parameters.insert<DenseMatrix<Real> *>(
      "Plane material matrix");
    equation_systems.parameters.insert<DenseMatrix<Real> *>(
      "Plate material matrix");
    equation_systems.parameters.set<DenseMatrix<Real> *>(
      "Plane material matrix") = &Dm;
    equation_systems.parameters.set<DenseMatrix<Real> *>(
      "Plate material matrix") = &Dp;
  }

  /**
   * Transformes the element from global to local space. It constructs a matrix
   * storing the local positions of the element's nodes and the transformation
   * matrix itself. A matrix storing the first partial derivatives is also
   * created. The area of the element is stored as well.
   * @param elem pointer to the element, that is transformed into local space
   * (in-param)
   * @param transUV reference, out-param, stores positions of the element's
   * nodes in local coordinate system
   * @param trafo reference, out-param, the transformation matrix from global to
   * local co-sys
   * @param dphi reference, out-param, stores the first partial derivatives of
   * the element in local space
   * @param area pointer, out-param, stores the area of the element
   */
  void shellsolid::initElement(const Elem **elem,
                               DenseMatrix<Real> &transUV,
                               DenseMatrix<Real> &trafo,
                               DenseMatrix<Real> &dphi,
                               Real *area)
  {
    // needed to differently treat different element types:
    ElemType type = (*elem)->type();
    // temporarily stores pointers to element's nodes:
    const Node *ndi = NULL, *ndj = NULL;
    // 'Node' can also be used as mathematical vector. U,V,W act as vectors for
    // the local coordinate system:
    Node U, V, W;

    if (type == TRI3) // three-node triangular element Tri-3
      {
        // transform arbirtrary 3D triangle to xy-plane with node A at origin
        // (implicitly):
        ndi = (*elem)->node_ptr(0); // node A
        ndj = (*elem)->node_ptr(1); // node B
        U = (*ndj) - (*ndi);        // U = B-A
        ndj = (*elem)->node_ptr(2); // node C
        V = (*ndj) - (*ndi);        // V = C-A

        transUV.resize(3, 2);
        for (int i = 0; i < 3; i++)
          { // node A lies in local origin (per definition)
            transUV(i, 0) = U(i); // node B in global coordinates (triangle
                                  // translated s.t. A lies in origin)
            transUV(i, 1) = V(i); // node C in global coordinates ( -"- )
          }
        /* transUV [ b_x, c_x ]
         *         [ b_y, c_y ]
         *         [ b_z, c_z ]
         */
        W = U.cross(V);
        // area of triangle is half the length of the cross product of U and V
        *area = 0.5 * W.size();

        U = U.unit();   // local x-axis unit vector
        W = W.unit();   // local z-axis unit vector, normal to triangle
        V = W.cross(U); // local y-axis unit vector (cross product of 2
                        // normalized vectors is automatically normalized)
      }
    else if (type == QUAD4) // four-node quadrilateral element Quad-4
      {
        // transform planar 3D quadrilateral to xy-plane:
        Node nI, nJ, nK, nL;
        ndi = (*elem)->node_ptr(0);            // node A
        ndj = (*elem)->node_ptr(1);            // node B
        nI = (*ndi) + 0.5 * ((*ndj) - (*ndi)); // nI = midpoint on edge AB
        ndi = (*elem)->node_ptr(2);            // node C
        nJ = (*ndj) + 0.5 * ((*ndi) - (*ndj)); // nJ = midpoint on edge BC
        ndj = (*elem)->node_ptr(3);            // node D
        nK = (*ndi) + 0.5 * ((*ndj) - (*ndi)); // nK = midpoint on edge CD
        ndi = (*elem)->node_ptr(0);            // node A (again)
        nL = (*ndj) + 0.5 * ((*ndi) - (*ndj)); // nL = midpoint on edge DA
        ndj = (*elem)->node_ptr(2);

        transUV.resize(3, 4);
        for (int i = 0; i < 4; i++)
          {
            ndi = (*elem)->node_ptr(i);
            transUV(0, i) = (*ndi)(0); // coord x in global coordinates
            transUV(1, i) = (*ndi)(1); // coord y in global coordinates
            transUV(2, i) = (*ndi)(2); // coord z in global coordinates
          }
        /* transUV [ a_x, b_x, c_x, d_x ]
         *         [ a_y, b_y, c_y, d_y ]
         *         [ a_z, b_z, c_z, d_z ]
         */

        U = nJ - nL;    // Vx
        U = U.unit();   // Vx normalized -> local x-axis unit vector
        W = nK - nI;    // Vr
        W = U.cross(W); // Vz = Vx x Vr
        W = W.unit();   // Vz normalized -> local z-axis unit vector
        V = W.cross(U); // Vy = Vz x Vx -> local y-axis unit vector
      }
    // at this point, the local axes unit vectors are stored in U, V and W
    trafo.resize(3, 3); // global to local transformation matrix
    for (int i = 0; i < 3; i++)
      {
        trafo(0, i) = U(i);
        trafo(1, i) = V(i);
        trafo(2, i) = W(i);
      }
    /* trafo [ u_x, u_y, u_z ]
     *       [ v_x, v_y, v_z ]
     *       [ w_x, w_y, w_z ]
     */

    // transform element's nodes to local coordinates and store results in
    // transUV:
    transUV.left_multiply(trafo);

    /*if (debug)
        {
        std::cout << "transUV:" << std::endl;
                transUV.print(std::cout);
        std::cout << std::endl << "trafo:" << std::endl;
                trafo.print(std::cout);
                std::cout << std::endl;
    }*/

    // calculate the partial derivaties; differently for the single element
    // types:
    if (type == TRI3)
      {
        dphi.resize(3, 2);
        dphi(0, 0) = -transUV(0, 0);                // x12 = x1-x2 = 0-x2 = -x2
        dphi(1, 0) = transUV(0, 1);                 // x31 = x3-x1 = x3-0 = x3
        dphi(2, 0) = transUV(0, 0) - transUV(0, 1); // x23 = x2-x3
        dphi(0, 1) =
          -transUV(1, 0); // y12 = y1-y2 = -y2 = 0 (stays zero, as node
                          // B and A lie on local x-axis and therefore)
        dphi(1, 1) = transUV(1, 1);                 // y31 = y3-y1 = y3-0 = y3
        dphi(2, 1) = transUV(1, 0) - transUV(1, 1); // y23 = y2-y3 = 0-y3 = -y3
      }
    else if (type == QUAD4)
      {
        dphi.resize(4, 2);
        dphi(0, 0) = transUV(0, 0) - transUV(0, 1); // x12 = x1-x2
        dphi(1, 0) = transUV(0, 1) - transUV(0, 2); // x23 = x2-x3
        dphi(2, 0) = transUV(0, 2) - transUV(0, 3); // x34 = x3-x4
        dphi(3, 0) = transUV(0, 3) - transUV(0, 0); // x41 = x4-x1

        dphi(0, 1) = transUV(1, 0) - transUV(1, 1); // y12 = y1-y2
        dphi(1, 1) = transUV(1, 1) - transUV(1, 2); // y23 = y2-y3
        dphi(2, 1) = transUV(1, 2) - transUV(1, 3); // y34 = y3-y4
        dphi(3, 1) = transUV(1, 3) - transUV(1, 0); // y41 = y4-y1

        *area = 0.0;
        // Gauss's area formula:
        // x_i*y_{i+1} - x_{i+1}*y_i = det((x_i,x_i+1),(y_i,y_i+1))
        for (int i = 0; i < 4; i++)
          *area += transUV(0, i) * transUV(1, (i + 1) % 4) -
                   transUV(0, (i + 1) % 4) * transUV(1, i);
        *area *= 0.5;
      }
  }

  /**
   * Constructs the stiffness matrix for the plane element component
   * @param type Type of the current element
   * @param transUV reference, in-param, local positions of element's nodes
   * @param dphi reference, in-param, partial derivatives of the element
   * @param area pointer, in-param, area of the element
   * @param Ke_m reference, out-param, stiffness matrix for plane element
   * component
   */
  void shellsolid::calcPlane(EquationSystems &es,
                             ElemType type,
                             DenseMatrix<Real> &transUV,
                             DenseMatrix<Real> &dphi,
                             Real *area,
                             DenseMatrix<Real> &Ke_m)
  { // Unpack the material matrices
    DenseMatrix<Real> &Dm =
      *(es.parameters.get<DenseMatrix<Real> *>("Plane material matrix"));
    DenseMatrix<Real> &Dp =
      *(es.parameters.get<DenseMatrix<Real> *>("Plate material matrix"));
    // Unpack the thickness
    Real thickness = es.parameters.get<Real>("Thickness");
    if (type == TRI3)
      {
        // construct strain-displacement matrix B
        DenseMatrix<Real> B_m(3, 6);
        B_m(0, 0) = dphi(2, 1);  //  y23
        B_m(0, 2) = dphi(1, 1);  //  y31
        B_m(0, 4) = dphi(0, 1);  //  y12
        B_m(1, 1) = -dphi(2, 0); // -x23
        B_m(1, 3) = -dphi(1, 0); // -x31
        B_m(1, 5) = -dphi(0, 0); // -x12
        B_m(2, 0) = -dphi(2, 0); // -x23
        B_m(2, 1) = dphi(2, 1);  //  y23
        B_m(2, 2) = -dphi(1, 0); // -x31
        B_m(2, 3) = dphi(1, 1);  //  y31
        B_m(2, 4) = -dphi(0, 0); // -x12
        B_m(2, 5) = dphi(0, 1);  //  y12
        B_m *= 1.0 / (2.0 * (*area));

        // Ke_m = t*A* B^T * Dm * B
        Ke_m = Dm;                         // Ke_m = 3x3
        Ke_m.right_multiply(B_m);          // Ke_m = 3x6
        Ke_m.left_multiply_transpose(B_m); // Ke_m = 6x6
        Ke_m *= thickness * (*area);       // considered thickness and area is
                                           // constant all over the element
      }
    else if (type == QUAD4)
      {
        // quadrature points definition:
        Real root = sqrt(1.0 / 3.0); // note: sqrt(3)/3 <=> sqrt(1/3)

        DenseMatrix<Real> B_m;        // strain-displacement-matrix
        DenseMatrix<Real> G(4, 8);    // temp matrix
        DenseMatrix<Real> J(2, 2);    // Jacobian
        DenseVector<Real> shapeQ4(4); // evaluation of shape functions
        DenseVector<Real> dhdr(4),
          dhds(
            4); // derivatives of shape functions wrt local coordinates r and s

        // we iterate over the 4 Gauss quadrature points (+- sqrt(1/3)) with
        // weight
        // 1
        Ke_m.resize(8, 8); // the resulting stiffness matrix
        for (int ii = 0; ii < 2; ii++)
          {
            Real r = pow(-1.0, ii) * root; // +/- root
            for (int jj = 0; jj < 2; jj++)
              {
                Real s = pow(-1.0, jj) * root; // +/- root

                shapeQ4(0) = 0.25 * (1 - r) * (1 - s);
                shapeQ4(1) = 0.25 * (1 + r) * (1 - s);
                shapeQ4(2) = 0.25 * (1 + r) * (1 + s);
                shapeQ4(3) = 0.25 * (1 - r) * (1 + s);

                dhdr(0) = -0.25 * (1 - s);
                dhdr(1) = 0.25 * (1 - s);
                dhdr(2) = 0.25 * (1 + s);
                dhdr(3) = -0.25 * (1 + s);

                dhds(0) = -0.25 * (1 - r);
                dhds(1) = -0.25 * (1 + r);
                dhds(2) = 0.25 * (1 + r);
                dhds(3) = 0.25 * (1 - r);

                J.resize(2, 2); // resizing automatically zero-s entries
                for (int i = 0; i < 4; i++)
                  {
                    J(0, 0) += dhdr(i) * transUV(0, i);
                    J(0, 1) += dhdr(i) * transUV(1, i);
                    J(1, 0) += dhds(i) * transUV(0, i);
                    J(1, 1) += dhds(i) * transUV(1, i);
                  }

                Real detjacob = J.det(); // Jacobian determinant

                B_m.resize(3, 4);
                B_m(0, 0) = J(1, 1);
                B_m(0, 1) = -J(0, 1);
                B_m(1, 2) = -J(1, 0);
                B_m(1, 3) = J(0, 0);
                B_m(2, 0) = -J(1, 0);
                B_m(2, 1) = J(0, 0);
                B_m(2, 2) = J(1, 1);
                B_m(2, 3) = -J(0, 1);
                B_m *= 1.0 / detjacob;

                for (int i = 0; i < 4; i++)
                  {
                    G(0, 2 * i) = dhdr(i);
                    G(1, 2 * i) = dhds(i);
                    G(2, 1 + 2 * i) = dhdr(i);
                    G(3, 1 + 2 * i) = dhds(i);
                  }

                // final step to get the strain-displacement-matrix B:
                B_m.right_multiply(G);

                // Ke_m = t * B^T * Dm * B * |J|
                DenseMatrix<Real> Ke_m_tmp;
                Ke_m_tmp = Dm;                         // Ke_m = 3x3
                Ke_m_tmp.left_multiply_transpose(B_m); // Ke_m = 8x8
                Ke_m_tmp.right_multiply(B_m);          // Ke_m = 3x8
                Ke_m_tmp *=
                  detjacob * thickness; // considered thickness and area is
                                        // constant all over the element

                Ke_m += Ke_m_tmp;
              }
          } // end of Gauss sampling for-loops
      }     // end of element type switch
  }

  /**
   * Constructs the stiffness matrix for the plate element component
   * @param type Type of the current element
   * @param dphi reference, in-param, partial derivatives of the element
   * @param area pointer, in-param, area of the element
   * @param Ke_p reference, out-param, stiffness matrix for plate element
   * component
   */
  void shellsolid::calcPlate(EquationSystems &es,
                             ElemType type,
                             DenseMatrix<Real> &dphi,
                             Real *area,
                             DenseMatrix<Real> &Ke_p)
  { // Unpack the material matrices
    DenseMatrix<Real> &Dm =
      *(es.parameters.get<DenseMatrix<Real> *>("Plane material matrix"));
    DenseMatrix<Real> &Dp =
      *(es.parameters.get<DenseMatrix<Real> *>("Plate material matrix"));
    // Unpack the thickness
    Real thickness = es.parameters.get<Real>("Thickness");
    DenseVector<Real> sidelen; // stores squared side lengths of the element

    if (type == TRI3)
      {
        std::vector<std::vector<double>> qps(3); // quadrature points
        for (unsigned int i = 0; i < qps.size(); i++)
          qps[i].resize(2);
        qps[0][0] = 1.0 / 6.0;
        qps[0][1] = 1.0 / 6.0;
        qps[1][0] = 2.0 / 3.0;
        qps[1][1] = 1.0 / 6.0;
        qps[2][0] = 1.0 / 6.0;
        qps[2][1] = 2.0 / 3.0;

        // squared side lengths
        sidelen.resize(3);
        sidelen(0) =
          pow(dphi(0, 0), 2.0) + pow(dphi(0, 1), 2.0); // side AB, x12^2 + y12^2
        sidelen(1) =
          pow(dphi(1, 0), 2.0) + pow(dphi(1, 1), 2.0); // side AC, x31^2 + y31^2
        sidelen(2) =
          pow(dphi(2, 0), 2.0) + pow(dphi(2, 1), 2.0); // side BC, x23^2 + y23^2

        Ke_p.resize(9, 9);
        for (unsigned int i = 0; i < qps.size(); i++)
          {
            DenseMatrix<Real>
              B; // strain-displacement-matrix
                 // construct B and evaluate it at the quadrature point
            evalBTri(es, sidelen, qps[i][0], qps[i][1], dphi, B);

            // construct auxiliary matrix Y:
            DenseMatrix<Real> Y(3, 3);
            Y(0, 0) = pow(dphi(2, 1), 2.0);
            Y(0, 1) = pow(dphi(1, 1), 2.0);
            Y(0, 2) = dphi(2, 1) * dphi(1, 1);
            Y(1, 0) = pow(dphi(2, 0), 2.0);
            Y(1, 1) = pow(dphi(1, 0), 2.0);
            Y(1, 2) = dphi(1, 0) * dphi(2, 0);
            Y(2, 0) = -2.0 * dphi(2, 0) * dphi(2, 1);
            Y(2, 1) = -2.0 * dphi(1, 0) * dphi(1, 0);
            Y(2, 2) = -dphi(2, 0) * dphi(1, 1) - dphi(1, 0) * dphi(2, 1);
            Y *= 1.0 / (4.0 * pow(*area, 2.0));

            DenseMatrix<Real> temp;
            temp = Dp;                       // temp = 3x3
            temp.right_multiply(Y);          // temp = 3x3
            temp.right_multiply(B);          // temp = 9x3
            temp.left_multiply_transpose(Y); // temp = 9x3
            temp.left_multiply_transpose(B); // temp = 9x9

            temp *= 1.0 / 6.0; // gauss-weight

            Ke_p += temp;
          }

        Ke_p *= 2.0 * (*area);
      }
    else if (type == QUAD4)
      {
        // squared side lengths:
        sidelen.resize(4);
        sidelen(0) =
          pow(dphi(0, 0), 2.0) + pow(dphi(0, 1), 2.0); // side AB, x12^2 + y12^2
        sidelen(1) =
          pow(dphi(1, 0), 2.0) + pow(dphi(1, 1), 2.0); // side BC, x23^2 + y23^2
        sidelen(2) =
          pow(dphi(2, 0), 2.0) + pow(dphi(2, 1), 2.0); // side CD, x34^2 + y34^2
        sidelen(3) =
          pow(dphi(3, 0), 2.0) + pow(dphi(3, 1), 2.0); // side DA, x41^2 + y41^2

        DenseMatrix<Real> Hcoeffs(5,
                                  4); // [ a_k, b_k, c_k, d_k, e_k ], k=5,6,7,8
        for (int i = 0; i < 4; i++)
          {
            Hcoeffs(0, i) = -dphi(i, 0) / sidelen(i);                    // a_k
            Hcoeffs(1, i) = 0.75 * dphi(i, 0) * dphi(i, 1) / sidelen(i); // b_k
            Hcoeffs(2, i) =
              (0.25 * pow(dphi(i, 0), 2.0) - 0.5 * pow(dphi(i, 1), 2.0)) /
              sidelen(i);                             // c_k
            Hcoeffs(3, i) = -dphi(i, 1) / sidelen(i); // d_k
            Hcoeffs(4, i) =
              (0.25 * pow(dphi(i, 1), 2.0) - 0.5 * pow(dphi(i, 0), 2.0)) /
              sidelen(i); // e_k
          }
        /*if (debug)
                {
            std::cout << "Hcoeffs:" << std::endl;
            Hcoeffs.print(std::cout);
            std::cout << std::endl;
        }*/

        Ke_p.resize(12, 12);

        // quadrature points definition:
        Real root = sqrt(1.0 / 3.0);
        DenseMatrix<Real> J(2, 2), Jinv(2, 2); // Jacobian and its inverse
        for (int ii = 0; ii < 2; ii++)
          {
            Real r = pow(-1.0, ii) * root; // +/- sqrt(1/3)
            for (int jj = 0; jj < 2; jj++)
              {
                Real s = pow(-1.0, jj) * root; // +/- sqrt(1/3)

                J(0, 0) =
                  (dphi(0, 0) + dphi(2, 0)) * s - dphi(0, 0) + dphi(2, 0);
                J(0, 1) =
                  (dphi(0, 1) + dphi(2, 1)) * s - dphi(0, 1) + dphi(2, 1);
                J(1, 0) =
                  (dphi(0, 0) + dphi(2, 0)) * r - dphi(1, 0) + dphi(3, 0);
                J(1, 1) =
                  (dphi(0, 1) + dphi(2, 1)) * r - dphi(1, 1) + dphi(3, 1);
                J *= 0.25;
                /*if (debug)
                                {
                    std::cout << "J:" << std::endl;
                    J.print(std::cout);
                    std::cout << std::endl;
                }*/
                Real det = J.det();
                // if (debug)
                //    std::cout << "|J| = " << det << std::endl;

                Jinv(0, 0) = J(1, 1);
                Jinv(0, 1) = -J(0, 1);
                Jinv(1, 0) = -J(1, 0);
                Jinv(1, 1) = J(0, 0);
                Jinv *= 1.0 / det;

                /*if (debug)
                                {
                    std::cout << "Jinv:" << std::endl;
                    Jinv.print(std::cout);
                    std::cout << std::endl;
                }*/
                DenseMatrix<Real> B;
                // construct strain-displacement-matrix B and evaluate it at the
                // current quadrature point:
                evalBQuad(es, Hcoeffs, r, s, Jinv, B);

                /*if (debug)
                                {
                    std::cout << "B:" << std::endl;
                    B.print(std::cout);
                    std::cout << std::endl;
                }*/
                DenseMatrix<Real> temp;
                temp = Dp;                       // temp = 3x3
                temp.left_multiply_transpose(B); // temp = 12x3
                temp.right_multiply(B);          // temp = 12x12
                temp *= det;

                Ke_p += temp;
              }
          } // end of quadrature point for-loops
      }     // end of element type switch
  }

  /*Construct the strain displacement matrix B for tri-3 plane element*/
  /*Required for stress calculation */
  void shellsolid::B_plane_tri(Real *area,
                               DenseMatrix<Real> &dphi,
                               DenseMatrix<Real> &out)
  {
    out.resize(3, 6);
    out(0, 0) = dphi(2, 1);  //  y23
    out(0, 2) = dphi(1, 1);  //  y31
    out(0, 4) = dphi(0, 1);  //  y12
    out(1, 1) = -dphi(2, 0); // -x23
    out(1, 3) = -dphi(1, 0); // -x31
    out(1, 5) = -dphi(0, 0); // -x12
    out(2, 0) = -dphi(2, 0); // -x23
    out(2, 1) = dphi(2, 1);  //  y23
    out(2, 2) = -dphi(1, 0); // -x31
    out(2, 3) = dphi(1, 1);  //  y31
    out(2, 4) = -dphi(0, 0); // -x12
    out(2, 5) = dphi(0, 1);  //  y12
    out *= 1.0 / (2.0 * (*area));
  }

  /*Construct the strain displacement matrix B for quad-4 plane element*/
  /*Required for stress calculation */
  void shellsolid::B_plane_quad(Real *area,
                                DenseMatrix<Real> &dphi,
                                DenseMatrix<Real> &transUV,
                                Real qp_x,
                                Real qp_y,
                                DenseMatrix<Real> &out)
  {
    DenseMatrix<Real> G(4, 8);    // temp matrix
    DenseMatrix<Real> J(2, 2);    // Jacobian
    DenseVector<Real> shapeQ4(4); // evaluation of shape functions
    DenseVector<Real> dhdr(4),
      dhds(4); // derivatives of shape functions wrt local coordinates r and s

    shapeQ4(0) = 0.25 * (1 - qp_x) * (1 - qp_y);
    shapeQ4(1) = 0.25 * (1 + qp_x) * (1 - qp_y);
    shapeQ4(2) = 0.25 * (1 + qp_x) * (1 + qp_y);
    shapeQ4(3) = 0.25 * (1 - qp_x) * (1 + qp_y);

    dhdr(0) = -0.25 * (1 - qp_y);
    dhdr(1) = 0.25 * (1 - qp_y);
    dhdr(2) = 0.25 * (1 + qp_y);
    dhdr(3) = -0.25 * (1 + qp_y);

    dhds(0) = -0.25 * (1 - qp_x);
    dhds(1) = -0.25 * (1 + qp_x);
    dhds(2) = 0.25 * (1 + qp_x);
    dhds(3) = 0.25 * (1 - qp_x);

    J.resize(2, 2); // resizing automatically zero-s entries
    for (int i = 0; i < 4; i++)
      {
        J(0, 0) += dhdr(i) * transUV(0, i);
        J(0, 1) += dhdr(i) * transUV(1, i);
        J(1, 0) += dhds(i) * transUV(0, i);
        J(1, 1) += dhds(i) * transUV(1, i);
      }

    Real detjacob = J.det(); // Jacobian determinant

    out.resize(3, 4);
    out(0, 0) = J(1, 1);
    out(0, 1) = -J(0, 1);
    out(1, 2) = -J(1, 0);
    out(1, 3) = J(0, 0);
    out(2, 0) = -J(1, 0);
    out(2, 1) = J(0, 0);
    out(2, 2) = J(1, 1);
    out(2, 3) = -J(0, 1);
    out *= 1.0 / detjacob;

    for (int i = 0; i < 4; i++)
      {
        G(0, 2 * i) = dhdr(i);
        G(1, 2 * i) = dhds(i);
        G(2, 1 + 2 * i) = dhdr(i);
        G(3, 1 + 2 * i) = dhds(i);
      }

    // final step to get the strain-displacement-matrix B:
    out.right_multiply(G);
  }

  /**
   * Constructs the strain-displacement-matrix B for the Tri-3 plate element at
   * the specified quadrature point.
   * @param C reference, in-param, vector containing the squared side lengths
   * @param L1 in-param, first triangle coordinate component
   * @param L2 in-param, second triangle coordinate component
   * @param dphi reference, in-param, partial derivatives of the element
   * @param out reference, out-param, the strain-displacement-matrix to be
   * constructed
   */
  void shellsolid::evalBTri(EquationSystems &es,
                            DenseVector<Real> &C,
                            Real L1,
                            Real L2,
                            DenseMatrix<Real> &dphi,
                            DenseMatrix<Real> &out)
  { // Unpack the material matrices
    DenseMatrix<Real> &Dm =
      *(es.parameters.get<DenseMatrix<Real> *>("Plane material matrix"));
    DenseMatrix<Real> &Dp =
      *(es.parameters.get<DenseMatrix<Real> *>("Plate material matrix"));
    // Unpack the thickness
    Real thickness = es.parameters.get<Real>("Thickness");
    out.resize(3, 9); // the future B

    Real mu1 = (C(0) - C(1)) / C(2);
    Real mu2 = (C(2) - C(0)) / C(1);
    Real mu3 = (C(1) - C(2)) / C(0);

    // some abbreviations to shorten the following terms
    Real L3 = 1 - L1 - L2;
    Real f13mu1 = 1 + 3 * mu1;
    Real f13mu2 = 1 + 3 * mu2;
    Real f13mu3 = 1 + 3 * mu3;
    Real f1m3mu3 = 1 - 3 * mu3;
    Real fm13mu2 = -1 + 3 * mu2;
    Real fm1m3mu3 = -1 - 3 * mu3;
    Real f1mmu1 = 1 - mu1;
    Real f1mmu2 = 1 - mu2;
    Real f1mmu3 = 1 - mu3;

    Real a = 3 * f1mmu3 * L1 - f13mu3 * L2 + f13mu3 * L3;
    Real b = 3 * f1mmu2 * L3 - f13mu2 * L1 + f13mu2 * L2;
    Real c = 3 * f1mmu1 * L2 - f13mu1 * L3 + f13mu1 * L1;

    // see page 38f of the thesis:
    // the following terms contains second order derivatives of the 9 shape
    // functions wrt the triangle coordinates L1 and L2
    out(0, 0) = 6 + L2 * (-4 - 2 * a) + 4 * f1m3mu3 * (L2 * L3 - L1 * L2) -
                12 * L1 + 2 * L2 * b + 8 * (L2 * L3 - L1 * L2);

    out(0, 1) = -dphi(1, 1) *
                  (-2 + 6 * L1 + 4 * L2 - L2 * b - 4 * L2 * L3 + 4 * L1 * L2) -
                dphi(0, 1) * (2 * L2 - L2 * a + L2 * L3 * 2 * f1m3mu3 -
                              L1 * L2 * 2 * f1m3mu3);

    out(0, 2) =
      dphi(1, 0) * (-2 + 6 * L1 + 4 * L2 - L2 * b - 4 * L2 * L3 + 4 * L1 * L2) +
      dphi(0, 0) *
        (2 * L2 - L2 * a + L2 * L3 * 2 * f1m3mu3 - L1 * L2 * 2 * f1m3mu3);

    out(0, 3) = -2 * L2 * c + 4 * f13mu1 * (L2 * L3 - L1 * L2) - 4 * L2 +
                2 * L2 * a + 4 * f1m3mu3 * (-L2 * L3 + L1 * L2);

    out(0, 4) =
      -dphi(0, 1) *
        (2 * L2 - L2 * a + L2 * L3 * 2 * f1m3mu3 - L1 * L2 * 2 * f1m3mu3) -
      dphi(2, 1) * (-L2 * c + L2 * L3 * 2 * f13mu1 - L1 * L2 * 2 * f13mu1);

    out(0, 5) =
      dphi(0, 0) *
        (2 * L2 - L2 * a + L2 * L3 * 2 * f1m3mu3 - L1 * L2 * 2 * f1m3mu3) +
      dphi(2, 0) * (-L2 * c + L2 * L3 * 2 * f13mu1 - L1 * L2 * 2 * f13mu1);

    out(0, 6) = -6 + 12 * L1 + 8 * L2 - 2 * L2 * b + 8 * (L1 * L2 - L2 * L3) +
                2 * L2 * c + 4 * f13mu1 * (L1 * L2 - L2 * L3);

    out(0, 7) =
      -dphi(2, 1) * (-L2 * c + L2 * L3 * 2 * f13mu1 - L1 * L2 * 2 * f13mu1) -
      dphi(1, 1) * (-4 + 6 * L1 + 4 * L2 - L2 * b - 4 * L2 * L3 + 4 * L1 * L2);

    out(0, 8) =
      dphi(2, 0) * (-L2 * c + L2 * L3 * 2 * f13mu1 - L1 * L2 * 2 * f13mu1) +
      dphi(1, 0) * (-4 + 6 * L1 + 4 * L2 - L2 * b - 4 * L2 * L3 + 4 * L1 * L2);

    out(1, 0) = -2 * L1 * a + 2 * L1 * L3 * 2 * fm1m3mu3 -
                2 * L1 * L2 * 2 * fm1m3mu3 - 4 * L1 + 2 * L1 * b -
                2 * L1 * L3 * 2 * fm13mu2 + 2 * L1 * L2 * 2 * fm13mu2;

    out(1, 1) = -dphi(1, 1) * (2 * L1 - 1 * L1 * b + 1 * L1 * L3 * 2 * fm13mu2 -
                               1 * L1 * L2 * 2 * fm13mu2) -
                dphi(0, 1) * (-1 * L1 * a + 1 * L1 * L3 * 2 * fm1m3mu3 -
                              1 * L1 * L2 * 2 * fm1m3mu3);

    out(1, 2) = dphi(1, 0) * (2 * L1 - 1 * L1 * b + 1 * L1 * L3 * 2 * fm13mu2 -
                              1 * L1 * L2 * 2 * fm13mu2) +
                dphi(0, 0) * (-1 * L1 * a + 1 * L1 * L3 * 2 * fm1m3mu3 -
                              1 * L1 * L2 * 2 * fm1m3mu3);

    out(1, 3) = 6 - 12 * L2 - 4 * L1 - 2 * L1 * c + 8 * L3 * L1 - 8 * L1 * L2 +
                2 * L1 * a - 2 * L1 * L3 * 2 * fm1m3mu3 +
                2 * L1 * L2 * 2 * fm1m3mu3;

    out(1, 4) = -dphi(0, 1) * (-1 * L1 * a + 1 * L1 * L3 * 2 * fm1m3mu3 -
                               1 * L1 * L2 * 2 * fm1m3mu3) -
                dphi(2, 1) * (-6 * L2 + 2 - 2 * L1 - 1 * L1 * c + 4 * L3 * L1 -
                              4 * L1 * L2);

    out(1, 5) = dphi(0, 0) * (-1 * L1 * a + 1 * L1 * L3 * 2 * fm1m3mu3 -
                              1 * L1 * L2 * 2 * fm1m3mu3) +
                dphi(2, 0) * (-6 * L2 + 2 - 2 * L1 - 1 * L1 * c + 4 * L3 * L1 -
                              4 * L1 * L2);

    out(1, 6) = -6 + 8 * L1 - 2 * L1 * b + 2 * L1 * L3 * 2 * fm13mu2 -
                2 * L1 * L2 * 2 * fm13mu2 + 12 * L2 + 2 * L1 * c - 8 * L3 * L1 +
                8 * L1 * L2;

    out(1, 7) = -dphi(2, 1) * (-6 * L2 + 4 - 2 * L1 - 1 * L1 * c + 4 * L3 * L1 -
                               4 * L1 * L2) -
                dphi(1, 1) * (2 * L1 - 1 * L1 * b + 1 * L1 * L3 * 2 * fm13mu2 -
                              1 * L1 * L2 * 2 * fm13mu2);

    out(1, 8) = dphi(2, 0) * (-6 * L2 + 4 - 2 * L1 - 1 * L1 * c + 4 * L3 * L1 -
                              4 * L1 * L2) +
                dphi(1, 0) * (2 * L1 - 1 * L1 * b + 1 * L1 * L3 * 2 * fm13mu2 -
                              1 * L1 * L2 * 2 * fm13mu2);

    out(2, 0) = 2 - 4 * L1 + L3 * a - L2 * a + L2 * L3 * 2 * fm1m3mu3 - L1 * a -
                L1 * L2 * 2 * fm1m3mu3 + L1 * L3 * 2 * f1m3mu3 -
                L1 * L2 * 2 * f1m3mu3 - 4 * L2 - L3 * b + L2 * b -
                L2 * L3 * 2 * fm13mu2 + L1 * b + L1 * L2 * 2 * fm13mu2 +
                4 * L3 * L1 - 4 * L1 * L2;

    out(2, 1) =
      -dphi(1, 1) * (-1 + 4 * L1 + 2 * L2 + 0.5 * L3 * b - 0.5 * L2 * b +
                     0.5 * L2 * L3 * 2 * fm13mu2 - 0.5 * L1 * b -
                     0.5 * L1 * L2 * 2 * fm13mu2 - 2 * L3 * L1 + 2 * L1 * L2) -
      dphi(0, 1) *
        (2 * L1 + 0.5 * L3 * a - 0.5 * L2 * a + 0.5 * L2 * L3 * 2 * fm1m3mu3 -
         0.5 * L1 * a - 0.5 * L1 * L2 * 2 * fm1m3mu3 +
         0.5 * L1 * L3 * 2 * f1m3mu3 - 0.5 * L1 * L2 * 2 * f1m3mu3);

    out(2, 2) =
      dphi(1, 0) * (-1 + 4 * L1 + 2 * L2 + 0.5 * L3 * b - 0.5 * L2 * b +
                    0.5 * L2 * L3 * 2 * fm13mu2 - 0.5 * L1 * b -
                    0.5 * L1 * L2 * 2 * fm13mu2 - 2 * L3 * L1 + 2 * L1 * L2) +
      dphi(0, 0) *
        (2 * L1 + 0.5 * L3 * a - 0.5 * L2 * a + 0.5 * L2 * L3 * 2 * fm1m3mu3 -
         0.5 * L1 * a - 0.5 * L1 * L2 * 2 * fm1m3mu3 +
         0.5 * L1 * L3 * 2 * f1m3mu3 - 0.5 * L1 * L2 * 2 * f1m3mu3);

    out(2, 3) = 2 - 4 * L2 + L3 * c - L2 * c + 4 * L2 * L3 - L1 * c -
                4 * L1 * L2 + L1 * L3 * 2 * f13mu1 - L1 * L2 * 2 * f13mu1 -
                4 * L1 - L3 * a + L2 * a + L1 * a - L2 * L3 * 2 * fm1m3mu3 +
                L1 * L2 * 2 * fm1m3mu3 - L1 * L3 * 2 * f1m3mu3 +
                L1 * L2 * 2 * f1m3mu3;

    out(2, 4) =
      -dphi(0, 1) *
        (2 * L1 + 0.5 * L3 * a - 0.5 * L2 * a + 0.5 * L2 * L3 * 2 * fm1m3mu3 -
         0.5 * L1 * a - 0.5 * L1 * L2 * 2 * fm1m3mu3 +
         0.5 * L1 * L3 * 2 * f1m3mu3 - 0.5 * L1 * L2 * 2 * f1m3mu3 - 1) -
      dphi(2, 1) *
        (-2 * L2 + 0.5 * L3 * c - 0.5 * L2 * c + 2 * L2 * L3 - 0.5 * L1 * c -
         2 * L1 * L2 + 0.5 * L1 * L3 * 2 * f13mu1 - 0.5 * L1 * L2 * 2 * f13mu1);

    out(2, 5) =
      dphi(0, 0) *
        (2 * L1 + 0.5 * L3 * a - 0.5 * L2 * a + 0.5 * L2 * L3 * 2 * fm1m3mu3 -
         0.5 * L1 * a - 0.5 * L1 * L2 * 2 * fm1m3mu3 +
         0.5 * L1 * L3 * 2 * f1m3mu3 - 0.5 * L1 * L2 * 2 * f1m3mu3 - 1) +
      dphi(2, 0) *
        (-2 * L2 + 0.5 * L3 * c - 0.5 * L2 * c + 2 * L2 * L3 - 0.5 * L1 * c -
         2 * L1 * L2 + 0.5 * L1 * L3 * 2 * f13mu1 - 0.5 * L1 * L2 * 2 * f13mu1);

    out(2, 6) = -4 + 8 * L1 + 8 * L2 + L3 * b - L2 * b + L2 * L3 * 2 * fm13mu2 -
                L1 * b - L1 * L2 * 2 * fm13mu2 - 4 * L3 * L1 + 8 * L1 * L2 -
                L3 * c + L2 * c - 4 * L2 * L3 + L1 * c - L1 * L3 * 2 * f13mu1 +
                L1 * L2 * 2 * f13mu1;

    out(2, 7) =
      -dphi(2, 1) * (-2 * L2 + 0.5 * L3 * c - 0.5 * L2 * c + 2 * L2 * L3 -
                     0.5 * L1 * c - 2 * L1 * L2 + 0.5 * L1 * L3 * 2 * f13mu1 -
                     0.5 * L1 * L2 * 2 * f13mu1 + 1) -
      dphi(1, 1) * (-2 + 4 * L1 + 2 * L2 + 0.5 * L3 * b - 0.5 * L2 * b +
                    0.5 * L2 * L3 * 2 * fm13mu2 - 0.5 * L1 * b -
                    0.5 * L1 * L2 * 2 * fm13mu2 - 2 * L3 * L1 + 2 * L1 * L2);

    out(2, 8) =
      dphi(2, 0) * (-2 * L2 + 0.5 * L3 * c - 0.5 * L2 * c + 2 * L2 * L3 -
                    0.5 * L1 * c - 2 * L1 * L2 + 0.5 * L1 * L3 * 2 * f13mu1 -
                    0.5 * L1 * L2 * 2 * f13mu1 + 1) +
      dphi(1, 0) * (-2 + 4 * L1 + 2 * L2 + 0.5 * L3 * b - 0.5 * L2 * b +
                    0.5 * L2 * L3 * 2 * fm13mu2 - 0.5 * L1 * b -
                    0.5 * L1 * L2 * 2 * fm13mu2 - 2 * L3 * L1 + 2 * L1 * L2);
    // the last row of the matrix must be multipled by 2 (this way, the upper
    // terms gets a bit shorter...)
    for (int i = 0; i < 9; i++)
      out(2, i) *= 2.0;
  }

  /**
   * Constructs the strain-displacement-matrix B for the Quad-4 plate element at
   * the specified quadrature point.
   * @param Hcoeffs reference, in-param, matrix containing coefficients
   * @param xi in-param, first local coordinate component
   * @param eta in-param, second local coordinate component
   * @param Jinv reference, in-param, inverse Jacobian matrix
   * @param out reference, out-param, the strain-displacement-matrix to be
   * constructed
   */
  void shellsolid::evalBQuad(EquationSystems &es,
                             DenseMatrix<Real> &Hcoeffs,
                             Real xi,
                             Real eta,
                             DenseMatrix<Real> &Jinv,
                             DenseMatrix<Real> &out)
  { // Unpack the material matrices
    DenseMatrix<Real> &Dm =
      *(es.parameters.get<DenseMatrix<Real> *>("Plane material matrix"));
    DenseMatrix<Real> &Dp =
      *(es.parameters.get<DenseMatrix<Real> *>("Plate material matrix"));
    // Unpack the thickness
    Real thickness = es.parameters.get<Real>("Thickness");
    out.resize(3, 12); // the future B

    // first order derivatives of the shape functions evaluated wrt to xi and
    // eta
    DenseVector<Real> N_xi(8), N_eta(8);
    N_xi(0) = 0.25 * (2.0 * xi + eta) * (1.0 - eta);
    N_xi(1) = 0.25 * (2.0 * xi - eta) * (1.0 - eta);
    N_xi(2) = 0.25 * (2.0 * xi + eta) * (1.0 + eta);
    N_xi(3) = 0.25 * (2.0 * xi - eta) * (1.0 + eta);
    N_xi(4) = -xi * (1.0 - eta);
    N_xi(5) = 0.5 * (1.0 - pow(eta, 2.0));
    N_xi(6) = -xi * (1.0 + eta);
    N_xi(7) = -0.5 * (1.0 - pow(eta, 2.0));

    N_eta(0) = 0.25 * (2.0 * eta + xi) * (1.0 - xi);
    N_eta(1) = 0.25 * (2.0 * eta - xi) * (1.0 + xi);
    N_eta(2) = 0.25 * (2.0 * eta + xi) * (1.0 + xi);
    N_eta(3) = 0.25 * (2.0 * eta - xi) * (1.0 - xi);
    N_eta(4) = -0.5 * (1.0 - pow(xi, 2.0));
    N_eta(5) = -eta * (1.0 + xi);
    N_eta(6) = 0.5 * (1.0 - pow(xi, 2.0));
    N_eta(7) = -eta * (1.0 - xi);

    // to make the code more readable, the indices gets replaced by letters
    int a = 0, b = 1, c = 2, d = 3, e = 4;
    int i5 = 0, i6 = 1, i7 = 2, i8 = 3;

    // see page 43ff. of the thesis
    DenseVector<Real> Hx_xi(12), Hy_xi(12), Hx_eta(12), Hy_eta(12);
    Hx_xi(0) = 1.5 * (Hcoeffs(a, i5) * N_xi(4) - Hcoeffs(a, i8) * N_xi(7));
    Hx_xi(1) = Hcoeffs(b, i5) * N_xi(4) + Hcoeffs(b, i8) * N_xi(7);
    Hx_xi(2) = N_xi(0) - Hcoeffs(c, i5) * N_xi(4) - Hcoeffs(c, i8) * N_xi(7);
    Hx_xi(3) = 1.5 * (Hcoeffs(a, i6) * N_xi(5) - Hcoeffs(a, i5) * N_xi(4));
    Hx_xi(4) = Hcoeffs(b, i6) * N_xi(5) + Hcoeffs(b, i5) * N_xi(4);
    Hx_xi(5) = N_xi(1) - Hcoeffs(c, i6) * N_xi(5) - Hcoeffs(c, i5) * N_xi(4);
    Hx_xi(6) = 1.5 * (Hcoeffs(a, i7) * N_xi(6) - Hcoeffs(a, i6) * N_xi(5));
    Hx_xi(7) = Hcoeffs(b, i7) * N_xi(6) + Hcoeffs(b, i6) * N_xi(5);
    Hx_xi(8) = N_xi(2) - Hcoeffs(c, i7) * N_xi(6) - Hcoeffs(c, i6) * N_xi(5);
    Hx_xi(9) = 1.5 * (Hcoeffs(a, i8) * N_xi(7) - Hcoeffs(a, i7) * N_xi(6));
    Hx_xi(10) = Hcoeffs(b, i8) * N_xi(7) + Hcoeffs(b, i7) * N_xi(6);
    Hx_xi(11) = N_xi(3) - Hcoeffs(c, i8) * N_xi(7) - Hcoeffs(c, i7) * N_xi(6);

    Hy_xi(0) = 1.5 * (Hcoeffs(d, i5) * N_xi(4) - Hcoeffs(d, i8) * N_xi(7));
    Hy_xi(1) = -N_xi(0) + Hcoeffs(e, i5) * N_xi(4) + Hcoeffs(e, i8) * N_xi(7);
    Hy_xi(2) = -Hx_xi(1);
    Hy_xi(3) = 1.5 * (Hcoeffs(d, i6) * N_xi(5) - Hcoeffs(d, i5) * N_xi(4));
    Hy_xi(4) = -N_xi(1) + Hcoeffs(e, i6) * N_xi(5) + Hcoeffs(e, i5) * N_xi(4);
    Hy_xi(5) = -Hx_xi(4);
    Hy_xi(6) = 1.5 * (Hcoeffs(d, i7) * N_xi(6) - Hcoeffs(d, i6) * N_xi(5));
    Hy_xi(7) = -N_xi(2) + Hcoeffs(e, i7) * N_xi(6) + Hcoeffs(e, i6) * N_xi(5);
    Hy_xi(8) = -Hx_xi(7);
    Hy_xi(9) = 1.5 * (Hcoeffs(d, i8) * N_xi(7) - Hcoeffs(d, i7) * N_xi(6));
    Hy_xi(10) = -N_xi(3) + Hcoeffs(e, i8) * N_xi(7) + Hcoeffs(e, i7) * N_xi(6);
    Hy_xi(11) = -Hx_xi(10);

    Hx_eta(0) = 1.5 * (Hcoeffs(a, i5) * N_eta(4) - Hcoeffs(a, i8) * N_eta(7));
    Hx_eta(1) = Hcoeffs(b, i5) * N_eta(4) + Hcoeffs(b, i8) * N_eta(7);
    Hx_eta(2) =
      N_eta(0) - Hcoeffs(c, i5) * N_eta(4) - Hcoeffs(c, i8) * N_eta(7);
    Hx_eta(3) = 1.5 * (Hcoeffs(a, i6) * N_eta(5) - Hcoeffs(a, i5) * N_eta(4));
    Hx_eta(4) = Hcoeffs(b, i6) * N_eta(5) + Hcoeffs(b, i5) * N_eta(4);
    Hx_eta(5) =
      N_eta(1) - Hcoeffs(c, i6) * N_eta(5) - Hcoeffs(c, i5) * N_eta(4);
    Hx_eta(6) = 1.5 * (Hcoeffs(a, i7) * N_eta(6) - Hcoeffs(a, i6) * N_eta(5));
    Hx_eta(7) = Hcoeffs(b, i7) * N_eta(6) + Hcoeffs(b, i6) * N_eta(5);
    Hx_eta(8) =
      N_eta(2) - Hcoeffs(c, i7) * N_eta(6) - Hcoeffs(c, i6) * N_eta(5);
    Hx_eta(9) = 1.5 * (Hcoeffs(a, i8) * N_eta(7) - Hcoeffs(a, i7) * N_eta(6));
    Hx_eta(10) = Hcoeffs(b, i8) * N_eta(7) + Hcoeffs(b, i7) * N_eta(6);
    Hx_eta(11) =
      N_eta(3) - Hcoeffs(c, i8) * N_eta(7) - Hcoeffs(c, i7) * N_eta(6);

    Hy_eta(0) = 1.5 * (Hcoeffs(d, i5) * N_eta(4) - Hcoeffs(d, i8) * N_eta(7));
    Hy_eta(1) =
      -N_eta(0) + Hcoeffs(e, i5) * N_eta(4) + Hcoeffs(e, i8) * N_eta(7);
    Hy_eta(2) = -Hx_eta(1);
    Hy_eta(3) = 1.5 * (Hcoeffs(d, i6) * N_eta(5) - Hcoeffs(d, i5) * N_eta(4));
    Hy_eta(4) =
      -N_eta(1) + Hcoeffs(e, i6) * N_eta(5) + Hcoeffs(e, i5) * N_eta(4);
    Hy_eta(5) = -Hx_eta(4);
    Hy_eta(6) = 1.5 * (Hcoeffs(d, i7) * N_eta(6) - Hcoeffs(d, i6) * N_eta(5));
    Hy_eta(7) =
      -N_eta(2) + Hcoeffs(e, i7) * N_eta(6) + Hcoeffs(e, i6) * N_eta(5);
    Hy_eta(8) = -Hx_eta(7);
    Hy_eta(9) = 1.5 * (Hcoeffs(d, i8) * N_eta(7) - Hcoeffs(d, i7) * N_eta(6));
    Hy_eta(10) =
      -N_eta(3) + Hcoeffs(e, i8) * N_eta(7) + Hcoeffs(e, i7) * N_eta(6);
    Hy_eta(11) = -Hx_eta(10);

    // the final construction process of B:
    for (int i = 0; i < 12; i++)
      {
        out(0, i) = Jinv(0, 0) * Hx_xi(i) + Jinv(0, 1) * Hx_eta(i);
        out(1, i) = Jinv(1, 0) * Hy_xi(i) + Jinv(1, 1) * Hy_eta(i);
        out(2, i) = Jinv(0, 0) * Hy_xi(i) + Jinv(0, 1) * Hy_eta(i) +
                    Jinv(1, 0) * Hx_xi(i) + Jinv(1, 1) * Hx_eta(i);
      }
  }

  /**
   * Superimposes the plane and plate stiffness matrices to the shell stiffness
   * matrix.
   * @param type Type of the current element
   * @param Ke_m reference, in-param, stiffness matrix of plane element
   * component
   * @param Ke_p reference, in-param, stiffness matrix of plate element
   * component
   * @param K_out reference, out-param, resulting stiffness matrix of shell
   * element
   */
  void shellsolid::constructStiffnessMatrix(EquationSystems &es,
                                            ElemType type,
                                            DenseMatrix<Real> &Ke_m,
                                            DenseMatrix<Real> &Ke_p,
                                            DenseMatrix<Real> &K_out)
  { // Unpack the material matrices
    DenseMatrix<Real> &Dm =
      *(es.parameters.get<DenseMatrix<Real> *>("Plane material matrix"));
    DenseMatrix<Real> &Dp =
      *(es.parameters.get<DenseMatrix<Real> *>("Plate material matrix"));
    // Unpack the thickness
    Real thickness = es.parameters.get<Real>("Thickness");
    int nodes = 4; // predefine QUAD-4 element
    if (type == TRI3)
      nodes = 3;
    else if (type == QUAD4)
      nodes = 4;

    // size of stiffness matrix depends on number of nodes the element has
    K_out.resize(6 * nodes, 6 * nodes);

    // copy values from the (nodes X nodes) sub-matrices into shell element
    // matrix:
    for (int i = 0; i < nodes; i++)
      {
        for (int j = 0; j < nodes; j++)
          {
            // submatrix K_ij [6x6]
            K_out(6 * i, 6 * j) = Ke_m(2 * i, 2 * j);                 // uu
            K_out(6 * i, 6 * j + 1) = Ke_m(2 * i, 2 * j + 1);         // uv
            K_out(6 * i + 1, 6 * j) = Ke_m(2 * i + 1, 2 * j);         // vu
            K_out(6 * i + 1, 6 * j + 1) = Ke_m(2 * i + 1, 2 * j + 1); // vv
            K_out(2 + 6 * i, 2 + 6 * j) = Ke_p(3 * i, 3 * j);         // ww
            K_out(2 + 6 * i, 2 + 6 * j + 1) = Ke_p(3 * i, 3 * j + 1); // wx
            K_out(2 + 6 * i, 2 + 6 * j + 2) = Ke_p(3 * i, 3 * j + 2); // wy
            K_out(2 + 6 * i + 1, 2 + 6 * j) = Ke_p(3 * i + 1, 3 * j); // xw
            K_out(2 + 6 * i + 1, 2 + 6 * j + 1) =
              Ke_p(3 * i + 1, 3 * j + 1); // xx
            K_out(2 + 6 * i + 1, 2 + 6 * j + 2) =
              Ke_p(3 * i + 1, 3 * j + 2);                             // xy
            K_out(2 + 6 * i + 2, 2 + 6 * j) = Ke_p(3 * i + 2, 3 * j); // yw
            K_out(2 + 6 * i + 2, 2 + 6 * j + 1) =
              Ke_p(3 * i + 2, 3 * j + 1); // yx
            K_out(2 + 6 * i + 2, 2 + 6 * j + 2) =
              Ke_p(3 * i + 2, 3 * j + 2); // yy
          }
      }

    // for the sixths d.o.f. we need the an approximate value
    // for this value the maximum value of the diagonal entries
    // of each sub-matrix must be known:
    Real max_value;
    for (int zi = 0; zi < nodes; zi++)
      {
        for (int zj = 0; zj < nodes; zj++)
          {
            // search for max value in plane-matrix
            max_value = Ke_m(2 * zi, 2 * zj); // begin with uu value
            max_value =
              std::max(max_value, Ke_m(2 * zi + 1, 2 * zj + 1)); // test for vv
            // search for max value in plate-matrix
            max_value =
              std::max(max_value, Ke_p(3 * zi, 3 * zj)); // test for ww
            max_value = std::max(
              max_value, Ke_p(3 * zi + 1, 3 * zj + 1)); // test for t_x t_x
            max_value = std::max(
              max_value, Ke_p(3 * zi + 2, 3 * zj + 2)); // test for t_y t_y
            // take max from both and divide it by 1000
            max_value /= 1000.0;
            // set it at corresponding place
            K_out(5 + 6 * zi, 5 + 6 * zj) = max_value;
          }
      }
  }

  /**
   * Transforms the local shell stiffness matrix back to global space
   * @param type Type of the current element
   * @param trafo reference, in-param, transformation matrix
   * @param Ke_inout reference, inout-param, modifies the shell stiffness matrix
   * with global space entries
   */
  void shellsolid::localToGlobalTrafo(EquationSystems &es,
                                      ElemType type,
                                      DenseMatrix<Real> &trafo,
                                      DenseMatrix<Real> &Ke_inout)
  { // Unpack the material matrices
    DenseMatrix<Real> &Dm =
      *(es.parameters.get<DenseMatrix<Real> *>("Plane material matrix"));
    DenseMatrix<Real> &Dp =
      *(es.parameters.get<DenseMatrix<Real> *>("Plate material matrix"));
    // Unpack the thickness
    Real thickness = es.parameters.get<Real>("Thickness");
    int nodes = 4; // predefine with QUAD-4 element
    if (type == TRI3)
      nodes = 3;
    else if (type == QUAD4)
      nodes = 4;

    DenseMatrix<Real> KeSub(6, 6); // one of the (nodes X nodes) sub-matrices
    DenseMatrix<Real> KeNew(
      6 * nodes, 6 * nodes); // the global version of the stiffness matrix
    DenseMatrix<Real> TSub(
      6, 6); // the transformation matrix for a stiffness sub-matrix
             // copy trafo two times into TSub (cf. comment below)
    for (int k = 0; k < 2; k++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          TSub(3 * k + i, 3 * k + j) = trafo(i, j);
    /* TSub: [ux, vx, wx,  0,  0,  0]
     *       [uy, vy, wy,  0,  0,  0]
     *       [uz, vz, wz,  0,  0,  0]
     *       [0 ,  0,  0, ux, vx, wx]
     *       [0 ,  0,  0, uy, vy, wy]
     *       [0 ,  0,  0, uz, vz, wz] */

    for (int i = 0; i < nodes; i++)
      {
        for (int j = 0; j < nodes; j++)
          {
            // copy values into temporary sub-matrix for correct format to
            // transformation
            for (int k = 0; k < 6; k++)
              for (int l = 0; l < 6; l++)
                KeSub(k, l) = Ke_inout(i * 6 + k, j * 6 + l);

            // the actual transformation step
            KeSub.right_multiply(TSub);
            KeSub.left_multiply_transpose(TSub);

            // copy transformed values into new global stiffness matrix
            for (int k = 0; k < 6; k++)
              for (int l = 0; l < 6; l++)
                KeNew(i * 6 + k, j * 6 + l) = KeSub(k, l);
          }
      }

    // bring stiffness matrix into right format for libMesh equation system
    // handling
    for (int alpha = 0; alpha < 6; alpha++)
      for (int beta = 0; beta < 6; beta++)
        for (int i = 0; i < nodes; i++)
          for (int j = 0; j < nodes; j++)
            Ke_inout(nodes * alpha + i, nodes * beta + j) =
              KeNew(6 * i + alpha, 6 * j + beta);
  }

  /**
   * Constructs the RHS for the given element by contributing the correct force
   * values for the so far unprocessed nodes.
   * @param elem pointer to the element whose right-hand side (RHS) is assembled
   * @param Fe reference, out-param, stores the RHS values
   * @param processedNodes pointer, inout-param, contains the already processed
   * nodes of the mesh to prevent double contribution of a force to the RHS
   */
  void shellsolid::contribRHS(EquationSystems &es,
                              const Elem **elem,
                              DenseVector<Real> &Fe,
                              std::unordered_set<unsigned int> *processedNodes)
  { // Unpack the material matrices
    DenseMatrix<Real> &Dm =
      *(es.parameters.get<DenseMatrix<Real> *>("Plane material matrix"));
    DenseMatrix<Real> &Dp =
      *(es.parameters.get<DenseMatrix<Real> *>("Plate material matrix"));
    // Unpack the thickness
    Real thickness = es.parameters.get<Real>("Thickness");
    // Unpack debug option
    bool debug = es.parameters.get<bool>("Debug");

    std::vector<DenseVector<Real>> &forces =
      (*es.parameters.get<std::vector<DenseVector<Real>> *>("Forcing terms"));
    unsigned int nsides =
      (*elem)->n_sides();  // 'sides' in libMesh equals 'nodes' for an element
    Fe.resize(6 * nsides); // prepare the element's RHS

    // go through all nodes (sides) of the element
    for (unsigned int side = 0; side < nsides; side++)
      {
        const Node *node =
          (*elem)->node_ptr(side);   // pointer to the current node
        dof_id_type id = node->id(); // ID of the current mesh node
        // do not process nodes that are owned by another process (for parallel
        // mode)
        if (node->processor_id() != global_processor_id())
          continue;

        // do not process already processed nodes
        if (processedNodes->find(id) == processedNodes->end())
          {
            // we process it now, so mark the node as processed for the future
            processedNodes->insert(id);

            for (int i = 0; i < 6; i++)
              Fe(side + nsides * i) = forces[id](i);
            // u_i   (0, 1,   ..., n-1)
            // v_i   (n, n+1, ...,2n-1)
            // w_i   (2n,2n+1,...,3n-1)
            // M_x_i (3n,3n+1,...,4n-1)
            // M_y_i (4n,4n+1,...,5n-1)
            // M_z_i (5n,5n+1,...,6n-1)
            if (debug)
              std::cout << "Forces for node " << id << ": Fx= " << Fe(side)
                        << ", Fy= " << Fe(side + nsides)
                        << ", Fz= " << Fe(side + nsides * 2)
                        << ", Mx= " << Fe(side + nsides * 3)
                        << ", My= " << Fe(side + nsides * 4)
                        << ", Mz= " << Fe(side + nsides * 5) << std::endl;
          }
      }
  }

  /**
   * Called by libMesh just before solving the system. The system matrix and RHS
   * will be assembled.
   * @param es reference, in-param, EquationSystem to get access to structures
   * like the mesh and the system
   * @param system_name reference, in-param The name of the system to assemble
   * the system matrix and RHS
   */
  void shellsolid::assemble_elasticity(EquationSystems &es,
                                       const std::string &system_name)
  {
    // Unpack the material matrices
    DenseMatrix<Real> &Dm =
      *(es.parameters.get<DenseMatrix<Real> *>("Plane material matrix"));
    DenseMatrix<Real> &Dp =
      *(es.parameters.get<DenseMatrix<Real> *>("Plate material matrix"));
    // Unpack the thickness
    Real thickness = es.parameters.get<Real>("Thickness");
    // Unpack the forcing term
    std::vector<DenseVector<Real>> &forces =
      *(es.parameters.get<std::vector<DenseVector<Real>> *>("Forcing terms"));

    // only allow the call for the Elasticity system
    libmesh_assert_equal_to(system_name, "Elasticity");

    // get a reference to the mesh
    const MeshBase &mesh = es.get_mesh();

    // get a reference to the linear implicit system
    LinearImplicitSystem &system =
      es.get_system<LinearImplicitSystem>("Elasticity");

    // A reference to the DofMap object for this system.
    // The DofMap object handles the index translation from node and element
    // numbers to degree of freedom numbers
    const DofMap &dof_map = system.get_dof_map();

    // stiffness matrices: Ke   for shell element,
    //                     Ke_m for plane element component ('m' like membrane),
    //                     Ke_p for plate element component ('p' like plate):
    DenseMatrix<Number> Ke, Ke_m, Ke_p;
    // RHS / force-momentum-vector:
    DenseVector<Number> Fe;

    // indices (positions) of node's variables in system matrix and RHS:
    std::vector<dof_id_type> dof_indices;

    DenseMatrix<Real>
      trafo; // global to local coordinate system transformation matrix
    DenseMatrix<Real>
      transUV; // stores the transformed positions of the element's nodes
    DenseMatrix<Real>
      dphi;          // contains the first partial derivatives of the element
    Real area = 0.0; // the area of the element

    // every node must contribute only once to the RHS. Since a node can be
    // shared by many elements 'processedNodes' keeps track of already used
    // nodes and prevent further processing of those
    std::unordered_set<dof_id_type> processedNodes;
    // we only need only as many nodes as the process has in his mesh partition
    processedNodes.reserve(mesh.n_local_nodes());

    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();
    // go through for all local elements in the mesh:
    for (; el != end_el; ++el)
      {
        const Elem *elem = *el;

        // get the local to global DOF-mappings for this element
        dof_map.dof_indices(elem, dof_indices);

        // get the type of the element
        ElemType type = elem->type();

        // transform the element from global to local space
        initElement(&elem, transUV, trafo, dphi, &area);

        // construct the stiffness matrices for the plane and plate element
        // component
        calcPlane(es, type, transUV, dphi, &area, Ke_m);
        calcPlate(es, type, dphi, &area, Ke_p);

        // superimpose both stiffness matrices to the shell element matrix
        constructStiffnessMatrix(es, type, Ke_m, Ke_p, Ke);

        // transform the shell stiffness matrix from local back to global space
        localToGlobalTrafo(es, type, trafo, Ke);

        // construct the right-hand side for the element
        contribRHS(es, &elem, Fe, &processedNodes);

        // constrain the matrix and RHS based on the defined boundary conditions
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        // add the element's matrix and RHS to the overall system matrix
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
      }
  }

  /* Stress Calculation*/
  void shellsolid::stress_calculation()
  {
    const unsigned int dim = mesh.mesh_dimension();

    unsigned int displacement_vars[6];
    displacement_vars[0] = system.variable_number("u");
    displacement_vars[1] = system.variable_number("v");
    displacement_vars[2] = system.variable_number("w");
    displacement_vars[3] = system.variable_number("tx");
    displacement_vars[4] = system.variable_number("ty");
    displacement_vars[5] = system.variable_number("tz");

    const DofMap &dof_map = system.get_dof_map();
    std::vector<std::vector<dof_id_type>> dof_indices_var(
      system.n_vars()); // Vector stores global dof dof_indices_var[var
                        // no.][node no.]

    unsigned int sigma_vars[3][3];
    sigma_vars[0][0] = stress_system.variable_number("sigma_xx");
    sigma_vars[0][1] = stress_system.variable_number("sigma_xy");
    sigma_vars[0][2] = stress_system.variable_number("sigma_xz");
    sigma_vars[1][0] = stress_system.variable_number("sigma_yx");
    sigma_vars[1][1] = stress_system.variable_number("sigma_yy");
    sigma_vars[1][2] = stress_system.variable_number("sigma_yz");
    sigma_vars[2][0] = stress_system.variable_number("sigma_zx");
    sigma_vars[2][1] = stress_system.variable_number("sigma_zy");
    sigma_vars[2][2] = stress_system.variable_number("sigma_zz");

    unsigned int sigma_vars_b[3][3];
    sigma_vars_b[0][0] = stress_system_b.variable_number("sigma_b_xx");
    sigma_vars_b[0][1] = stress_system_b.variable_number("sigma_b_xy");
    sigma_vars_b[0][2] = stress_system_b.variable_number("sigma_b_xz");
    sigma_vars_b[1][0] = stress_system_b.variable_number("sigma_b_yx");
    sigma_vars_b[1][1] = stress_system_b.variable_number("sigma_b_yy");
    sigma_vars_b[1][2] = stress_system_b.variable_number("sigma_b_yz");
    sigma_vars_b[2][0] = stress_system_b.variable_number("sigma_b_zx");
    sigma_vars_b[2][1] = stress_system_b.variable_number("sigma_b_zy");
    sigma_vars_b[2][2] = stress_system_b.variable_number("sigma_b_zz");

    const DofMap &stress_dof_map = stress_system.get_dof_map();
    std::vector<dof_id_type> stress_dof_indices_var;
    const DofMap &stress_b_dof_map = stress_system_b.get_dof_map();
    std::vector<dof_id_type> stress_b_dof_indices_var;
    DenseMatrix<Real> B_m; // Strain displacement matrix
    DenseMatrix<Real> B_b; // Strain displacement matrix bending
    DenseMatrix<Real>
      trafo; // global to local coordinate system transformation matrix
    DenseMatrix<Real> trafo_2; // tranformation matrix 6x6
    DenseMatrix<Real>
      transUV; // stores the transformed positions of the element's nodes
    DenseMatrix<Real>
      dphi;          // contains the first partial derivatives of the element
    Real area = 0.0; // the area of the element

    // plane element stress calculation
    DenseVector<Real>
      elem_dis_local; // element displacement at eact node (u1, v1, w1, tx1,
                      // tx2, tx3 ...) in local co-ordinates
    DenseVector<Real> elem_dis_local_uv; // only u and v componenets required
                                         // for stress calculation
    DenseVector<Real> elem_stress; // stress vector sigma_x, sigma_y, tau_xy
                                   // (local co-ordinates)
    DenseMatrix<Real> temp;
    DenseMatrix<Real> elem_stress_tensor; // plane element stress tensor

    // plate bending stress calculation
    DenseVector<Real> elem_dis_local_w; // only w, tx, ty componenets required
                                        // for bending stress calculation
    DenseVector<Real> elem_stress_b;    // bending stress vector M_x/I, M_y/I,
                                        // M_xy/I (local co-ordinates)
    DenseMatrix<Real> elem_stress_tensor_b; // bending stress tensor
    DenseVector<Real> sidelen; // stores squared side lengths of the element

    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();
    for (; el != end_el; ++el)
      {
        const Elem *elem = *el;
        ElemType type = elem->type();

        for (unsigned int var = 0; var < 6; var++)
          {
            dof_map.dof_indices(
              elem, dof_indices_var[var], displacement_vars[var]);
          }
        initElement(&elem, transUV, trafo, dphi, &area);
        // creating 6x6 transformation matrix
        trafo_2.resize(6, 6);
        for (int k = 0; k < 2; k++)
          {
            for (int i = 0; i < 3; i++)
              {
                for (int j = 0; j < 3; j++)
                  {
                    trafo_2(3 * k + i, 3 * k + j) = trafo(i, j);
                  }
              }
          }

        if (type == TRI3)
          {
            int nodes = 3;
            elem_stress.resize(3);
            temp.resize(3, 3);
            B_m.resize(3, 6);
            elem_dis_local.resize(6 * nodes);
            elem_dis_local_uv.resize(2 * nodes);
            elem_dis_local_w.resize(3 * nodes);
            elem_stress_tensor.resize(3, 3);
            elem_stress_tensor_b.resize(3, 3);

            B_plane_tri(&area, dphi, B_m);

            // transforming solution from global to local co-ordinates
            for (int i = 0; i < nodes; i++)
              {
                for (int j = 0; j < 6; j++)
                  {
                    for (int l = 0; l < 6; l++)
                      {
                        elem_dis_local(j + 6 * i) +=
                          trafo_2(j, l) *
                          system.current_solution(dof_indices_var[l][i]);
                      }
                  }
              }

            for (int i = 0; i < nodes; i++)
              {
                // selecting only u and v components required for plane stress
                // calculation
                elem_dis_local_uv(2 * i) = elem_dis_local(6 * i);
                elem_dis_local_uv(2 * i + 1) = elem_dis_local(6 * i + 1);
                // selecting only w, tx, ty components required for bending
                // stress calculation
                elem_dis_local_w(3 * i) = elem_dis_local(6 * i + 2);
                elem_dis_local_w(3 * i + 1) = elem_dis_local(6 * i + 3);
                elem_dis_local_w(3 * i + 2) = elem_dis_local(6 * i + 4);
              }

            temp = Dm;
            temp.right_multiply(B_m);

            // plane stress calculation
            for (unsigned int i = 0; i < 3; i++)
              {
                for (unsigned int j = 0; j < 6; j++)
                  {
                    elem_stress(i) += temp(i, j) * elem_dis_local_uv(j);
                  }
              }

            elem_stress_tensor(0, 0) = elem_stress(0);
            elem_stress_tensor(0, 1) = elem_stress(2);
            elem_stress_tensor(1, 0) = elem_stress(2);
            elem_stress_tensor(1, 1) = elem_stress(1);

            // Transforming stress tensor to global co-ordinate system
            elem_stress_tensor.right_multiply(trafo);
            elem_stress_tensor.left_multiply_transpose(trafo);

            // quadrature points
            std::vector<std::vector<double>> qps(3);
            for (unsigned int i = 0; i < qps.size(); i++)
              {
                qps[i].resize(2);
              }
            qps[0][0] = 1.0 / 6.0;
            qps[0][1] = 1.0 / 6.0;
            qps[1][0] = 2.0 / 3.0;
            qps[1][1] = 1.0 / 6.0;
            qps[2][0] = 1.0 / 6.0;
            qps[2][1] = 2.0 / 3.0;

            // squared side lengths (required for calculation of strain
            // diaplacement matrix bending)
            sidelen.resize(3);
            sidelen(0) = pow(dphi(0, 0), 2.0) +
                         pow(dphi(0, 1), 2.0); // side AB, x12^2 + y12^2
            sidelen(1) = pow(dphi(1, 0), 2.0) +
                         pow(dphi(1, 1), 2.0); // side AC, x31^2 + y31^2
            sidelen(2) = pow(dphi(2, 0), 2.0) +
                         pow(dphi(2, 1), 2.0); // side BC, x23^2 + y23^2

            for (unsigned int q = 0; q < qps.size(); q++)
              {
                elem_stress_b.resize(3);
                B_b.resize(3, 9);
                // construct B and evaluate it at the quadrature point
                evalBTri(
                  equation_systems, sidelen, qps[q][0], qps[q][1], dphi, B_b);
                temp = Dp;
                temp.right_multiply(B_b);
                // bending stress calculation
                for (unsigned int i = 0; i < 3; i++)
                  {
                    for (unsigned int j = 0; j < 3 * nodes; j++)
                      {
                        elem_stress_b(i) += temp(i, j) * elem_dis_local_w(j);
                      }
                  }

                elem_stress_tensor_b(0, 0) =
                  elem_stress_tensor_b(0, 0) + elem_stress_b(0);
                elem_stress_tensor_b(0, 1) =
                  elem_stress_tensor_b(0, 1) + elem_stress_b(2);
                elem_stress_tensor_b(1, 0) =
                  elem_stress_tensor_b(1, 0) + elem_stress_b(2);
                elem_stress_tensor_b(1, 1) =
                  elem_stress_tensor_b(1, 1) + elem_stress_b(1);
              }
            elem_stress_tensor_b.scale(
              1.0 / 3.0); // Averaging stress at quadrature points
            // Moment of inertia hard coded for 10x10 square plate
            double m_inertia = 10.0 * pow(thickness, 3.0) / 12.0;
            elem_stress_tensor_b.scale(1.0 / m_inertia);
            // Transforming stress tensor to global co-ordinate system
            elem_stress_tensor_b.right_multiply(trafo);
            elem_stress_tensor_b.left_multiply_transpose(trafo);
          }

        if (type == QUAD4)
          {
            int nodes = 4;
            elem_dis_local.resize(6 * nodes);
            elem_dis_local_uv.resize(2 * nodes);
            elem_dis_local_w.resize(3 * nodes);
            elem_stress_tensor.resize(3, 3);
            elem_stress_tensor_b.resize(3, 3);

            for (int i = 0; i < nodes; i++)
              {
                for (int j = 0; j < 6; j++)
                  {
                    for (int k = 0; k < 6; k++)
                      {
                        elem_dis_local(j + 6 * i) +=
                          trafo_2(j, k) *
                          system.current_solution(dof_indices_var[k][i]);
                      }
                  }
              }

            for (int i = 0; i < nodes; i++)
              {
                elem_dis_local_uv(2 * i) = elem_dis_local(6 * i);
                elem_dis_local_uv(2 * i + 1) = elem_dis_local(6 * i + 1);
                elem_dis_local_w(3 * i) = elem_dis_local(6 * i + 2);
                elem_dis_local_w(3 * i + 1) = elem_dis_local(6 * i + 3);
                elem_dis_local_w(3 * i + 2) = elem_dis_local(6 * i + 4);
              }

            // squared side lengths:
            sidelen.resize(4);
            sidelen(0) = pow(dphi(0, 0), 2.0) +
                         pow(dphi(0, 1), 2.0); // side AB, x12^2 + y12^2
            sidelen(1) = pow(dphi(1, 0), 2.0) +
                         pow(dphi(1, 1), 2.0); // side BC, x23^2 + y23^2
            sidelen(2) = pow(dphi(2, 0), 2.0) +
                         pow(dphi(2, 1), 2.0); // side CD, x34^2 + y34^2
            sidelen(3) = pow(dphi(3, 0), 2.0) +
                         pow(dphi(3, 1), 2.0); // side DA, x41^2 + y41^2

            DenseMatrix<Real> Hcoeffs(
              5, 4); // [ a_k, b_k, c_k, d_k, e_k ], k=5,6,7,8
            for (int i = 0; i < 4; i++)
              {
                Hcoeffs(0, i) = -dphi(i, 0) / sidelen(i); // a_k
                Hcoeffs(1, i) =
                  0.75 * dphi(i, 0) * dphi(i, 1) / sidelen(i); // b_k
                Hcoeffs(2, i) =
                  (0.25 * pow(dphi(i, 0), 2.0) - 0.5 * pow(dphi(i, 1), 2.0)) /
                  sidelen(i);                             // c_k
                Hcoeffs(3, i) = -dphi(i, 1) / sidelen(i); // d_k
                Hcoeffs(4, i) =
                  (0.25 * pow(dphi(i, 1), 2.0) - 0.5 * pow(dphi(i, 0), 2.0)) /
                  sidelen(i); // e_k
              }

            DenseMatrix<Real> J(2, 2), Jinv(2, 2); // Jacobian and its inverse
            Real root = sqrt(1.0 / 3.0);           // Quadrature point
            for (int ii = 0; ii < 2; ii++)
              {
                Real r = pow(-1.0, ii) * root; // +/- root
                for (int jj = 0; jj < 2; jj++)
                  {
                    Real s = pow(-1.0, jj) * root; // +/- root

                    elem_stress.resize(3);
                    elem_stress_b.resize(3);
                    temp.resize(3, 3);
                    B_m.resize(3, 8);
                    B_b.resize(3, 12);

                    J(0, 0) =
                      (dphi(0, 0) + dphi(2, 0)) * s - dphi(0, 0) + dphi(2, 0);
                    J(0, 1) =
                      (dphi(0, 1) + dphi(2, 1)) * s - dphi(0, 1) + dphi(2, 1);
                    J(1, 0) =
                      (dphi(0, 0) + dphi(2, 0)) * r - dphi(1, 0) + dphi(3, 0);
                    J(1, 1) =
                      (dphi(0, 1) + dphi(2, 1)) * r - dphi(1, 1) + dphi(3, 1);
                    J *= 0.25;
                    Real det = J.det();

                    Jinv(0, 0) = J(1, 1);
                    Jinv(0, 1) = -J(0, 1);
                    Jinv(1, 0) = -J(1, 0);
                    Jinv(1, 1) = J(0, 0);
                    Jinv *= 1.0 / det;

                    // construct strain-displacement-matrix B and evaluate it at
                    // the current quadrature point:
                    B_plane_quad(&area, dphi, transUV, r, s, B_m);
                    evalBQuad(equation_systems, Hcoeffs, r, s, Jinv, B_b);

                    temp = Dm;
                    temp.right_multiply(B_m);

                    for (unsigned int i = 0; i < 3; i++)
                      {
                        for (unsigned int j = 0; j < 2 * nodes; j++)
                          {
                            elem_stress(i) += temp(i, j) * elem_dis_local_uv(j);
                          }
                      }

                    elem_stress_tensor(0, 0) =
                      elem_stress_tensor(0, 0) + elem_stress(0);
                    elem_stress_tensor(0, 1) =
                      elem_stress_tensor(0, 1) + elem_stress(2);
                    elem_stress_tensor(1, 0) =
                      elem_stress_tensor(1, 0) + elem_stress(2);
                    elem_stress_tensor(1, 1) =
                      elem_stress_tensor(1, 1) + elem_stress(1);

                    temp = Dp;
                    temp.right_multiply(B_b);

                    for (unsigned int i = 0; i < 3; i++)
                      {
                        for (unsigned int j = 0; j < 3 * nodes; j++)
                          {
                            elem_stress_b(i) +=
                              temp(i, j) * elem_dis_local_w(j);
                          }
                      }

                    elem_stress_tensor_b(0, 0) =
                      elem_stress_tensor_b(0, 0) + elem_stress_b(0);
                    elem_stress_tensor_b(0, 1) =
                      elem_stress_tensor_b(0, 1) + elem_stress_b(2);
                    elem_stress_tensor_b(1, 0) =
                      elem_stress_tensor_b(1, 0) + elem_stress_b(2);
                    elem_stress_tensor_b(1, 1) =
                      elem_stress_tensor_b(1, 1) + elem_stress_b(1);
                  }
              }
            // Averaging stress at quadrature points
            elem_stress_tensor.scale(1.0 / 4.0);
            elem_stress_tensor_b.scale(1.0 / 4.0);
            // Moment of inertia hard coded for 10x10 square plate
            double m_inertia = 10.0 * pow(thickness, 3.0) / 12.0;
            elem_stress_tensor_b.scale(1.0 / m_inertia);
            // Transforming stress tensor to global co-ordinate system
            elem_stress_tensor.right_multiply(trafo);
            elem_stress_tensor.left_multiply_transpose(trafo);
            elem_stress_tensor_b.right_multiply(trafo);
            elem_stress_tensor_b.left_multiply_transpose(trafo);
          }

        // saving element stress in the system
        for (unsigned int i = 0; i < 3; i++)
          {
            for (unsigned int j = 0; j < 3; j++)
              {
                stress_dof_map.dof_indices(
                  elem, stress_dof_indices_var, sigma_vars[i][j]);
                dof_id_type dof_index = stress_dof_indices_var[0];

                if ((stress_system.solution->first_local_index() <=
                     dof_index) &&
                    (dof_index < stress_system.solution->last_local_index()))
                  {
                    stress_system.solution->set(dof_index,
                                                elem_stress_tensor(i, j));
                  }

                stress_b_dof_map.dof_indices(
                  elem, stress_b_dof_indices_var, sigma_vars_b[i][j]);
                dof_id_type dof_index_b = stress_b_dof_indices_var[0];

                if ((stress_system_b.solution->first_local_index() <=
                     dof_index_b) &&
                    (dof_index_b <
                     stress_system_b.solution->last_local_index()))
                  {
                    stress_system_b.solution->set(dof_index_b,
                                                  elem_stress_tensor_b(i, j));
                  }
              }
          }
      }

    stress_system.solution->close();
    stress_system.update();
    stress_system_b.solution->close();
    stress_system_b.update();
  }

  /**
   * creates output file of the displaced mesh with all solution data
   * @param mesh reference, in-param, needed as parameter for mesh export object
   * @param es reference, in-param, from this system collection the solution is
   * written to the file
   */
  void shellsolid::writeOutput()
  {
    // write the solution to file
    std::ostringstream file_name;
    file_name << "solid.e";

    ExodusII_IO(mesh).write_equation_systems(file_name.str(), equation_systems);
  }
} // namespace ShellSolid
