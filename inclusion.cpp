
#include "inclusion_solver.hpp"
#include "mfem.hpp"
#include <fstream>
#include <iostream>

#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include <Omega_h_adapt.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_timer.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_vector.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_curve_coarsen.hpp>

#include "miniapps/meshing/mesh-optimizer.hpp"
namespace oh = Omega_h;

namespace { // anonymous namespace

template <oh::Int dim>
static void set_target_metric(oh::Mesh* mesh, oh::Int scale, ParOmegaMesh
  *pOmesh, oh::Real error_des2) {
  auto coords = mesh->coords();
  auto target_metrics_w = oh::Write<oh::Real>
    (mesh->nverts() * oh::symm_ncomps(dim));
  pOmesh->ProjectFieldElementtoVertex (mesh, "zz_error");

  auto hd_hc = oh::Write<oh::Real>(mesh->nverts());

  auto error_c = mesh->get_array<oh::Real> (0, "zz_error");

  auto ev2v = mesh->get_adj(1,0).ab2b;
  auto length_edg_w = oh::Write<oh::Real> (mesh->nedges());
  auto f1 = OMEGA_H_LAMBDA(oh::LO e) {
    auto v0 = ev2v[e*2 + 0];
    auto v1 = ev2v[e*2 + 1];
    auto p0 = oh::get_vector<dim>(coords, v0);
    auto p1 = oh::get_vector<dim>(coords, v1);
    oh::Real dist = 0.0;
    for (oh::Int i = 0; i < dim; ++i) {
      dist += (p1[i] - p0[i])*(p1[i] - p0[i]);
    }
    dist = std::pow(dist, 0.5);
    length_edg_w[e] = dist;
  };
  oh::parallel_for(mesh->nedges(), f1);
  mesh->add_tag(oh::EDGE, "length_parent", 1, oh::Reals(length_edg_w));
  oh::ProjectFieldtoVertex (mesh, "length_parent", 1);
  auto length_c = mesh->get_array<oh::Real> (0, "length_parent");

  auto v_class_id = mesh->get_array<oh::LO> (0, "class_id");

  auto f = OMEGA_H_LAMBDA(oh::LO v) {
    auto h = oh::Vector<dim>();
    auto vtxError = error_c[v];
    for (oh::Int i = 0; i < dim; ++i) {
      h[i] = std::pow((error_des2/vtxError), 0.5)*length_c[v];
      if ((v_class_id[v] == 190) ||  (v_class_id[v] == 186)) {
        h[i] = length_c[v];
        //printf("inclusion vert %d\n", v);
      }
      hd_hc[v] = h[i]/length_c[v];
    }
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  oh::parallel_for(mesh->nverts(), f);
  mesh->set_tag(oh::VERT, "target_metric", oh::Reals(target_metrics_w));
  mesh->add_tag(oh::VERT, "hd_hc", 1, oh::Reals(hd_hc));
}

template <oh::Int dim>
void run_case(oh::Mesh* mesh, char const* vtk_path, oh::Int scale,
              const oh::Int myid, ParOmegaMesh *pOmesh, const oh::Real error_des2) {
  printf("in run case\n");
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = get_implied_metrics(mesh);
  mesh->add_tag(oh::VERT, "metric", oh::symm_ncomps(dim), implied_metrics);
  mesh->add_tag<oh::Real>(oh::VERT, "target_metric", oh::symm_ncomps(dim));
  set_target_metric<dim>(mesh, scale, pOmesh, error_des2);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->ask_lengths();
  mesh->ask_qualities();
  oh::vtk::FullWriter writer;
  if (vtk_path) {
    writer = oh::vtk::FullWriter(vtk_path, mesh);
    writer.write();
  }
  auto opts = oh::AdaptOpts(mesh);
  opts.should_swap = false;
  opts.should_filter_invalids = -1;
  //opts.should_coarsen = false;
  //opts.should_coarsen_slivers = false;
  opts.xfer_opts.type_map["zz_error"] = OMEGA_H_POINTWISE;
  opts.min_quality_allowed = 0.01;
  //opts.min_quality_desired = 0.1;
  //opts.max_length_allowed = 4.0*opts.max_length_desired;
  oh::Now t0 = oh::now();
 
  for (int itr=0; itr<1; ++itr) {
  //while (approach_metric(mesh, opts)) {
    printf("approach metric %d\n",approach_metric(mesh, opts));
    adapt(mesh, opts);
    if (vtk_path) {
      writer.write();
    }
  }
 
  //adapt(mesh, opts);
  oh::Now t1 = oh::now();
  if (!myid) std::cout << "total adapt time: " << (t1 - t0) << " seconds\n";

}

} // end anonymous namespace

// Permittivity Functions
Coefficient *
SetupPermittivityCoefficient(int max_attr, // max attribute in the mesh
    double eps1, double eps2, // Permittivity of each phase
    const Array<int>& p1, const Array<int>& p2); // list of regions for each phase

// Phi Boundary Condition
double phi_bc_uniform(const Vector &);

int main(int argc, char *argv[]) {
  MPI_Session mpi(argc, argv);
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Read Omega_h mesh
  auto lib = oh::Library();
  oh::Mesh o_mesh(&lib);
  //oh::binary::read ("../meshes/setup_1x1_3mil.osh", lib.world(), &o_mesh);
  oh::binary::read ("../meshes/setup_1x1_crv-coarse.osh", lib.world(), &o_mesh);
  
  if (o_mesh.is_curved() > 0) {
    printf("elevating to order 3\n");
    oh::calc_quad_ctrlPts_from_interpPts(&o_mesh);
    oh::elevate_curve_order_2to3(&o_mesh);
    o_mesh.add_tag<oh::Real>(0, "bezier_pts", o_mesh.dim(), o_mesh.coords());
    printf("checking validity of initial mesh\n");
    check_validity_all_tet(&o_mesh);
    auto cubic_curveVtk = oh::Mesh(o_mesh.comm()->library());
    cubic_curveVtk.set_comm(lib.world());
    build_cubic_curveVtk_3d(&o_mesh, &cubic_curveVtk, 5);
    std::string vtuPath = "/lore/joshia5/Meshes/curved/incl_ini_curveVtk.vtu";
    vtuPath += to_string(myid);
    oh::vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk, 2);
    auto cubic_wireframe = oh::Mesh(o_mesh.comm()->library());
    cubic_wireframe.set_comm(lib.world());
    build_cubic_wireframe_3d(&o_mesh, &cubic_wireframe, 5);
    vtuPath = "/lore/joshia5/Meshes/curved/incl_ini_wireframe.vtu";
    vtuPath += to_string(myid);
    oh::vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_wireframe, 1);
  }

  //Device device("cuda");
  //has very little effect, can be default for cuda build
  //device.Print();
  double error_des = 0.0;

  int max_iter = 2;

  // NOTE: the list containing model tags are model dependent
  double kappa = 2.; // relative permittivity of phase 2 wrt vacuum
  double epsilon1 = epsilon0_; // permittivity of substrate phase (1)
  double epsilon2 = kappa * epsilon0_; // permittivity of inclusion phase (2)

  const char *mesh_file = "";
  int order = o_mesh.get_max_order();
  bool visualization = false;

  Array<int> dbcs;
  Array<int> phase1; // list of model regions containing the phase 1 dielectric
  Array<int> phase2; // list of model regions containing the phase 2 dielectric


  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh",
      "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
      "Finite element order (polynomial degree).");
  args.AddOption(&phase1, "-s", "--substrate",
      "List of Model Regions Containing Phase 1 (Substrate).");
  args.AddOption(&phase2, "-i", "--inclusion",
      "List of Model Regions Containing Phase 2 (Inclusion).");
  args.AddOption(&kappa, "-k", "--kappa",
      "Relative permittivity phase 2.");
  args.AddOption(&dbcs, "-dbcs", "--dirichlet-bc-surf",
      "Dirichlet Boundary Condition Surfaces");
  args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
      "--no-visualization",
      "Enable or disable GLVis visualization.");
  args.Parse();
  if (!args.Good()) {
    if (mpi.Root()) {
      args.PrintUsage(cout);
    }
    return 1;
  }
  if (mpi.Root()) {
    args.PrintOptions(cout);
  }

  for (int Itr = 0; Itr < max_iter; Itr++)  {
    Mesh *mesh = new OmegaMesh (&o_mesh);
    ParMesh *pmesh = new ParMesh (MPI_COMM_WORLD, *mesh);
    //ParMesh *pmesh = new ParOmegaMesh (MPI_COMM_WORLD, &o_mesh);
    cout << "is NC" << pmesh->Nonconforming();

    int dim = pmesh->SpaceDimension();

    if (mpi.Root()) {
      cout << "Starting initialization." << endl;
    }

    // Create a coefficient describing the dielectric permittivity
    oh::Now t0 = oh::now();
    Coefficient * epsCoef =
      SetupPermittivityCoefficient(pmesh->attributes.Max(),
          epsilon1, epsilon2, phase1, phase2);

    // Boundary Conditions
    // Boundary conditions in this example are all around Dirichlet of the form
    // phi = -E.x, where E is a constant vector. This type of boundary conditions
    // mimics a uniform background electric field. Thus boundary conditions are
    // provided by the function phi_bc_uniform(x).

    // Create the Electrostatic solver
    InclusionSolver Inclusion(*pmesh, order, dbcs, *epsCoef, phi_bc_uniform);

    // Initialize GLVis visualization
    if (visualization) {
      Inclusion.InitializeGLVis();
    }

    if (mpi.Root()) { cout << "Initialization done." << endl; }

    // Display the current number of DoFs in each finite element space
    Inclusion.PrintSizes();

    // Assemble all forms
    Inclusion.Assemble();

    // Solve the system and compute any auxiliary fields
    Inclusion.Solve();

    oh::Now t1 = oh::now();
    if (!myid) std::cout << "total solve time: " << (t1 - t0) << " seconds\n";

    // Determine the current size of the linear system
    int prob_size = Inclusion.GetProblemSize();

    // Send the solution by socket to a GLVis server.
    if (visualization) {
      Inclusion.DisplayToGLVis();
    }

    Vector errors(pmesh->GetNE());
    Inclusion.GetErrorEstimates(errors);
    
    // adapt
    char Fname[128];
    sprintf(Fname,
        "inclusion_1x1_crv.vtk");
    char iter_str[8];
    sprintf(iter_str, "_%d", Itr);
    strcat(Fname, iter_str);
    puts(Fname);

    ParOmegaMesh* pOmesh = dynamic_cast<ParOmegaMesh*>(pmesh);
    pOmesh->ElementFieldMFEMtoOmegaH (&o_mesh, errors, dim, "zz_error");

    for (int i = 0; i < pmesh->GetNE(); i++) {
      //printf("elem %d , vol %1.10f, err %1.10f\n", i,pmesh->GetElementVolume(i),errors[i]); 
    }
    double local_max_err = errors.Max();
    double local_min_err = errors.Min();
    //printf("before adapt run case error max %1.10f error min %1.10f \n",
    //  local_max_err, local_min_err);
    double global_max_err;
    MPI_Allreduce(&local_max_err, &global_max_err, 1,
        MPI_DOUBLE, MPI_MAX, pmesh->GetComm());

    double global_min_err;
    MPI_Allreduce(&local_min_err, &global_min_err, 1,
        MPI_DOUBLE, MPI_MIN, pmesh->GetComm());

    double errorVol_tot_loc = 0.0;
    double error_tot_loc = 0.0;
    double vol_tot_loc = 0.0;
    for (int i = 0; i < pmesh->GetNE(); i++) {
      auto vol_loc = pmesh->GetElementVolume(i);
      vol_tot_loc += vol_loc;
      errorVol_tot_loc += errors[i]*vol_loc;
      error_tot_loc += errors[i];
    }
    double vol_tot;
    double errorVol_tot;
    double error_tot;
    MPI_Allreduce(&vol_tot_loc, &vol_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&errorVol_tot_loc, &errorVol_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&error_tot_loc, &error_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double error_bar = errorVol_tot/vol_tot;
    //const double frac = 40;
    const double frac = 0.9;
    if (Itr == 0) error_des = frac*error_bar;

    printf("Itr %d error tot %1.10f error_bar %1.15f #tet %d tot_vol %1.10f\n",
        Itr, error_tot, error_bar, mesh->GetNE(), vol_tot);
    //if (error_bar < error_des) break;

    stringstream ss;
    ss << "inclusion_iter_" << Itr;
    Inclusion.WriteToVtk(ss.str().c_str());

    if (global_max_err < error_des) {
      oh::vtk::write_parallel("before_break", &o_mesh);
      cout << "converged\n";
      break;
    }
    
    {
      std::string mesh_path = "1x1_crv_";
      mesh_path += std::to_string(3);
      mesh_path += ".mesh";
      ofstream mesh_ofs(mesh_path);
      mesh_ofs.precision(8);
      pmesh->Print(mesh_ofs);
      std::string gf_path = "eField_";
      gf_path += std::to_string(3);
      gf_path += ".gf";
      ofstream eField_ofs(gf_path);
      eField_ofs.precision(8);
      Inclusion.e_->Save(eField_ofs);
    }
    /*
    */

    if (o_mesh.is_curved() > 0) {
      GridFunction *nodes = pmesh->GetNodes();
      ofstream nodes_ofs("nodes.gf");
      nodes_ofs.precision(8);
      nodes->Save(nodes_ofs);

      auto cubic_curveVtk = oh::Mesh(o_mesh.comm()->library());
      cubic_curveVtk.set_comm(lib.world());
      build_cubic_curveVtk_3d(&o_mesh, &cubic_curveVtk, 5);
      std::string vtuPath = "/lore/joshia5/Meshes/curved/incl_bef_curveVtk.vtu";
      oh::vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk, 2);
      auto cubic_wireframe = oh::Mesh(o_mesh.comm()->library());
      cubic_wireframe.set_comm(lib.world());
      build_cubic_wireframe_3d(&o_mesh, &cubic_wireframe, 5);
      vtuPath = "/lore/joshia5/Meshes/curved/incl_bef_wireframe.vtu";
      oh::vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_wireframe, 1);
    }

    if (Itr+1 < max_iter) {
      run_case<3>(&o_mesh, Fname, Itr, myid, pOmesh, error_des);
    }
    else {
      fprintf(stderr, "calling TMOP\n");
      {
        //Set inputs as per the 3D untangling case
        int mesh_poly_deg     = 3;
        int metric_id = 313;
        int target_id         = 1;
        double solver_rtol    = 1e-5;
        int max_lin_iter      = 50;
        int quad_order        = 4;
        bool fdscheme         = true;
        int verbosity_level   = 1;
        //defaults
        int barrier_type = 0;
        int worst_case_type    = 0;
        int h_metric_id       = -1;
        bool pa               = false;
        bool normalization    = false;
        double lim_const      = 0.0;
        double adapt_lim_const   = 0.0;
        int adapt_eval        = 0;
        int n_h_iter          = 1;
        int n_hr_iter         = 5;
        int solver_art_type   = 0;
        int solver_iter       = 20;
        bool surface_fit_adapt = false;
        double surface_fit_threshold = -10;
        int lin_solver        = 2;
        int solver_type       = 0;
        bool hradaptivity     = false;
        bool move_bnd         = true;
        int combomet          = 0;
        bool exactaction      = false;

        // 11. Form the integrator that uses the chosen metric and target.
        double min_detJ = -0.1;
        TMOP_QualityMetric *metric = NULL;
        // T-metrics
        metric = new TMOP_Metric_313(min_detJ);

        TMOP_QualityMetric *h_metric = NULL;
        

        /*
        TMOP_WorstCaseUntangleOptimizer_Metric::BarrierType btype;
        btype = TMOP_WorstCaseUntangleOptimizer_Metric::BarrierType::None;
        TMOP_WorstCaseUntangleOptimizer_Metric::WorstCaseType wctype;
        wctype = TMOP_WorstCaseUntangleOptimizer_Metric::WorstCaseType::None;
        HessianCoefficient *adapt_coeff = NULL;
        HRHessianCoefficient *hr_adapt_coeff = NULL;
        */

        TMOP_QualityMetric *untangler_metric = NULL;

        TargetConstructor::TargetType target_t;
        TargetConstructor *target_c = NULL;
        H1_FECollection ind_fec(mesh_poly_deg, dim);
        FiniteElementSpace ind_fes(mesh, &ind_fec);
        FiniteElementSpace ind_fesv(mesh, &ind_fec, dim);
        GridFunction size(&ind_fes), aspr(&ind_fes), ori(&ind_fes);
        GridFunction aspr3d(&ind_fesv);

        const AssemblyLevel al =
          pa ? AssemblyLevel::PARTIAL : AssemblyLevel::LEGACY;

        target_t = TargetConstructor::IDEAL_SHAPE_UNIT_SIZE;

        if (target_c == NULL)
        {
          target_c = new TargetConstructor(target_t);
        }
        target_c->SetNodes(x0);
        TMOP_QualityMetric *metric_to_use = barrier_type > 0 || worst_case_type > 0
          ? untangler_metric
          : metric;
        TMOP_Integrator *tmop_integ = new TMOP_Integrator(metric_to_use, target_c,
            h_metric);

        // Finite differences for computations of derivatives.
        if (fdscheme)
        {
          MFEM_VERIFY(pa == false, "PA for finite differences is not implemented.");
          tmop_integ->EnableFiniteDifferences(x);
        }
        tmop_integ->SetExactActionFlag(exactaction);

        // Setup the quadrature rules for the TMOP integrator.
        //int quad_type         = 1;
        IntegrationRules *irules = NULL;
        irules = &IntRulesLo;
        tmop_integ->SetIntegrationRules(*irules, quad_order);

        if (dim == 3)
        {
          cout << "Tetrahedron quadrature points: "
            << irules->Get(Geometry::TETRAHEDRON, quad_order).GetNPoints()
            << "\nHexahedron quadrature points: "
            << irules->Get(Geometry::CUBE, quad_order).GetNPoints()
            << "\nPrism quadrature points: "
            << irules->Get(Geometry::PRISM, quad_order).GetNPoints() << endl;
        }

        // Limit the node movement.
        // The limiting distances can be given by a general function of space.
        FiniteElementSpace dist_fespace(mesh, fec); // scalar space
        GridFunction dist(&dist_fespace);
        dist = 1.0;
        ConstantCoefficient lim_coeff(lim_const);
        if (lim_const != 0.0) { tmop_integ->EnableLimiting(x0, dist, lim_coeff); }

        // Adaptive limiting.
        GridFunction adapt_lim_gf0(&ind_fes);
        ConstantCoefficient adapt_lim_coeff(adapt_lim_const);
        AdaptivityEvaluator *adapt_lim_eval = NULL;
        if (adapt_lim_const > 0.0)
        {
          MFEM_VERIFY(pa == false, "PA is not implemented for adaptive limiting");

          FunctionCoefficient adapt_lim_gf0_coeff(adapt_lim_fun);
          adapt_lim_gf0.ProjectCoefficient(adapt_lim_gf0_coeff);

          if (adapt_eval == 0) { adapt_lim_eval = new AdvectorCG(al); }
          else if (adapt_eval == 1)
          {
#ifdef MFEM_USE_GSLIB
            adapt_lim_eval = new InterpolatorFP;
#else
            MFEM_ABORT("MFEM is not built with GSLIB support!");
#endif
          }
          else { MFEM_ABORT("Bad interpolation option."); }

          tmop_integ->EnableAdaptiveLimiting(adapt_lim_gf0, adapt_lim_coeff,
              *adapt_lim_eval);
          if (visualization)
          {
            socketstream vis1;
            common::VisualizeField(vis1, "localhost", 19916, adapt_lim_gf0, "Zeta 0",
                300, 600, 300, 300);
          }
        }

        // Surface fitting.
        L2_FECollection mat_coll(0, dim);
        H1_FECollection surf_fit_fec(mesh_poly_deg, dim);
        FiniteElementSpace surf_fit_fes(mesh, &surf_fit_fec);
        FiniteElementSpace mat_fes(mesh, &mat_coll);
        GridFunction mat(&mat_fes);
        GridFunction surf_fit_mat_gf(&surf_fit_fes);
        GridFunction surf_fit_gf0(&surf_fit_fes);
        Array<bool> surf_fit_marker(surf_fit_gf0.Size());
        ConstantCoefficient surf_fit_coeff(surface_fit_const);
        AdaptivityEvaluator *adapt_surface = NULL;
        if (surface_fit_const > 0.0)
        {
          MFEM_VERIFY(hradaptivity == false,
              "Surface fitting with HR is not implemented yet.");
          MFEM_VERIFY(pa == false,
              "Surface fitting with PA is not implemented yet.");

          FunctionCoefficient ls_coeff(surface_level_set);
          surf_fit_gf0.ProjectCoefficient(ls_coeff);

          for (int i = 0; i < mesh->GetNE(); i++)
          {
            mat(i) = material_id(i, surf_fit_gf0);
            mesh->SetAttribute(i, static_cast<int>(mat(i) + 1));
          }

          GridFunctionCoefficient mat_coeff(&mat);
          surf_fit_mat_gf.ProjectDiscCoefficient(mat_coeff, GridFunction::ARITHMETIC);
          for (int j = 0; j < surf_fit_marker.Size(); j++)
          {
            if (surf_fit_mat_gf(j) > 0.1 && surf_fit_mat_gf(j) < 0.9)
            {
              surf_fit_marker[j] = true;
              surf_fit_mat_gf(j) = 1.0;
            }
            else
            {
              surf_fit_marker[j] = false;
              surf_fit_mat_gf(j) = 0.0;
            }
          }

          if (adapt_eval == 0) { adapt_surface = new AdvectorCG; }
          else if (adapt_eval == 1)
          {
#ifdef MFEM_USE_GSLIB
            adapt_surface = new InterpolatorFP;
#else
            MFEM_ABORT("MFEM is not built with GSLIB support!");
#endif
          }
          else { MFEM_ABORT("Bad interpolation option."); }

          tmop_integ->EnableSurfaceFitting(surf_fit_gf0, surf_fit_marker,
              surf_fit_coeff, *adapt_surface);
          if (visualization)
          {
            socketstream vis1, vis2, vis3;
            common::VisualizeField(vis1, "localhost", 19916, surf_fit_gf0, "Level Set 0",
                300, 600, 300, 300);
            common::VisualizeField(vis2, "localhost", 19916, mat, "Materials",
                600, 600, 300, 300);
            common::VisualizeField(vis3, "localhost", 19916, surf_fit_mat_gf,
                "Dofs to Move",
                900, 600, 300, 300);
          }
        }

        // Has to be after the enabling of the limiting / alignment, as it computes
        // normalization factors for these terms as well.
        if (normalization) { tmop_integ->EnableNormalization(x0); }

        // 12. Setup the final NonlinearForm (which defines the integral of interest,
        //     its first and second derivatives). Here we can use a combination of
        //     metrics, i.e., optimize the sum of two integrals, where both are
        //     scaled by used-defined space-dependent weights. Note that there are no
        //     command-line options for the weights and the type of the second
        //     metric; one should update those in the code.
        NonlinearForm a(fespace);
        if (pa) { a.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
        ConstantCoefficient *metric_coeff1 = NULL;
        TMOP_QualityMetric *metric2 = NULL;
        TargetConstructor *target_c2 = NULL;
        FunctionCoefficient metric_coeff2(weight_fun);

        // Explicit combination of metrics.
        if (combomet > 0)
        {
          // First metric.
          metric_coeff1 = new ConstantCoefficient(1.0);
          tmop_integ->SetCoefficient(*metric_coeff1);

          // Second metric.
          if (dim == 2) { metric2 = new TMOP_Metric_077; }
          else          { metric2 = new TMOP_Metric_315; }
          TMOP_Integrator *tmop_integ2 = NULL;
          if (combomet == 1)
          {
            target_c2 = new TargetConstructor(
                TargetConstructor::IDEAL_SHAPE_EQUAL_SIZE);
            target_c2->SetVolumeScale(0.01);
            target_c2->SetNodes(x0);
            tmop_integ2 = new TMOP_Integrator(metric2, target_c2, h_metric);
            tmop_integ2->SetCoefficient(metric_coeff2);
          }
          else { tmop_integ2 = new TMOP_Integrator(metric2, target_c, h_metric); }
          tmop_integ2->SetIntegrationRules(*irules, quad_order);
          if (fdscheme) { tmop_integ2->EnableFiniteDifferences(x); }
          tmop_integ2->SetExactActionFlag(exactaction);

          TMOPComboIntegrator *combo = new TMOPComboIntegrator;
          combo->AddTMOPIntegrator(tmop_integ);
          combo->AddTMOPIntegrator(tmop_integ2);
          if (normalization) { combo->EnableNormalization(x0); }
          if (lim_const != 0.0) { combo->EnableLimiting(x0, dist, lim_coeff); }

          a.AddDomainIntegrator(combo);
        }
        else
        {
          a.AddDomainIntegrator(tmop_integ);
        }

        if (pa) { a.Setup(); }

        // Compute the minimum det(J) of the starting mesh.
        min_detJ = infinity();
        const int NE = mesh->GetNE();
        for (int i = 0; i < NE; i++)
        {
          const IntegrationRule &ir =
            irules->Get(fespace->GetFE(i)->GetGeomType(), quad_order);
          ElementTransformation *transf = mesh->GetElementTransformation(i);
          for (int j = 0; j < ir.GetNPoints(); j++)
          {
            transf->SetIntPoint(&ir.IntPoint(j));
            min_detJ = min(min_detJ, transf->Jacobian().Det());
          }
        }
        cout << "Minimum det(J) of the original mesh is " << min_detJ << endl;

        if (min_detJ < 0.0 && barrier_type == 0
            && metric_id != 22 && metric_id != 211 && metric_id != 252
            && metric_id != 311 && metric_id != 313 && metric_id != 352)
        {
          MFEM_ABORT("The input mesh is inverted! Try an untangling metric.");
        }
        if (min_detJ < 0.0)
        {
          MFEM_VERIFY(target_t == TargetConstructor::IDEAL_SHAPE_UNIT_SIZE,
              "Untangling is supported only for ideal targets.");

          const DenseMatrix &Wideal =
            Geometries.GetGeomToPerfGeomJac(fespace->GetFE(0)->GetGeomType());
          min_detJ /= Wideal.Det();

          // Slightly below minJ0 to avoid div by 0.
          min_detJ -= 0.01 * h0.Min();
        }

        // For HR tests, the energy is normalized by the number of elements.
        const double init_energy = a.GetGridFunctionEnergy(x) /
          (hradaptivity ? mesh->GetNE() : 1);
        //double init_metric_energy = init_energy;
        if (lim_const > 0.0 || adapt_lim_const > 0.0 || surface_fit_const > 0.0)
        {
          lim_coeff.constant = 0.0;
          adapt_lim_coeff.constant = 0.0;
          surf_fit_coeff.constant   = 0.0;
          //init_metric_energy = a.GetGridFunctionEnergy(x) /
            (hradaptivity ? mesh->GetNE() : 1);
          lim_coeff.constant = lim_const;
          adapt_lim_coeff.constant = adapt_lim_const;
          surf_fit_coeff.constant   = surface_fit_const;
        }

        // Visualize the starting mesh and metric values.
        // Note that for combinations of metrics, this only shows the first metric.
        if (visualization)
        {
          char title[] = "Initial metric values";
          vis_tmop_metric_s(mesh_poly_deg, *metric, *target_c, *mesh, title, 0);
        }

        // 13. Fix all boundary nodes, or fix only a given component depending on the
        //     boundary attributes of the given mesh. Attributes 1/2/3 correspond to
        //     fixed x/y/z components of the node. Attribute 4 corresponds to an
        //     entirely fixed node. Other boundary attributes do not affect the node
        //     movement boundary conditions.
        if (move_bnd == false)
        {
          Array<int> ess_bdr(mesh->bdr_attributes.Max());
          ess_bdr = 1;
          a.SetEssentialBC(ess_bdr);
        }
        else
        {
          int n = 0;
          for (int i = 0; i < mesh->GetNBE(); i++)
          {
            const int nd = fespace->GetBE(i)->GetDof();
            const int attr = mesh->GetBdrElement(i)->GetAttribute();
            MFEM_VERIFY(!(dim == 2 && attr == 3),
                "Boundary attribute 3 must be used only for 3D meshes. "
                "Adjust the attributes (1/2/3/4 for fixed x/y/z/all "
                "components, rest for free nodes), or use -fix-bnd.");
            if (attr == 1 || attr == 2 || attr == 3) { n += nd; }
            if (attr == 4) { n += nd * dim; }
          }
          Array<int> ess_vdofs(n);
          n = 0;
          for (int i = 0; i < mesh->GetNBE(); i++)
          {
            const int nd = fespace->GetBE(i)->GetDof();
            const int attr = mesh->GetBdrElement(i)->GetAttribute();
            fespace->GetBdrElementVDofs(i, vdofs);
            if (attr == 1) // Fix x components.
            {
              for (int j = 0; j < nd; j++)
              { ess_vdofs[n++] = vdofs[j]; }
            }
            else if (attr == 2) // Fix y components.
            {
              for (int j = 0; j < nd; j++)
              { ess_vdofs[n++] = vdofs[j+nd]; }
            }
            else if (attr == 3) // Fix z components.
            {
              for (int j = 0; j < nd; j++)
              { ess_vdofs[n++] = vdofs[j+2*nd]; }
            }
            else if (attr == 4) // Fix all components.
            {
              for (int j = 0; j < vdofs.Size(); j++)
              { ess_vdofs[n++] = vdofs[j]; }
            }
          }
          a.SetEssentialVDofs(ess_vdofs);
        }

        // 14. As we use the Newton method to solve the resulting nonlinear system,
        //     here we setup the linear solver for the system's Jacobian.
        Solver *S = NULL, *S_prec = NULL;
        const double linsol_rtol = 1e-12;
        if (lin_solver == 0)
        {
          S = new DSmoother(1, 1.0, max_lin_iter);
        }
        else if (lin_solver == 1)
        {
          CGSolver *cg = new CGSolver;
          cg->SetMaxIter(max_lin_iter);
          cg->SetRelTol(linsol_rtol);
          cg->SetAbsTol(0.0);
          cg->SetPrintLevel(verbosity_level >= 2 ? 3 : -1);
          S = cg;
        }
        else
        {
          MINRESSolver *minres = new MINRESSolver;
          minres->SetMaxIter(max_lin_iter);
          minres->SetRelTol(linsol_rtol);
          minres->SetAbsTol(0.0);
          if (verbosity_level > 2) { minres->SetPrintLevel(1); }
          minres->SetPrintLevel(verbosity_level == 2 ? 3 : -1);
          if (lin_solver == 3 || lin_solver == 4)
          {
            if (pa)
            {
              MFEM_VERIFY(lin_solver != 4, "PA l1-Jacobi is not implemented");
              auto js = new OperatorJacobiSmoother;
              js->SetPositiveDiagonal(true);
              S_prec = js;
            }
            else
            {
              auto ds = new DSmoother((lin_solver == 3) ? 0 : 1, 1.0, 1);
              ds->SetPositiveDiagonal(true);
              S_prec = ds;
            }
            minres->SetPreconditioner(*S_prec);
          }
          S = minres;
        }

        // Perform the nonlinear optimization.
        const IntegrationRule &ir =
          irules->Get(fespace->GetFE(0)->GetGeomType(), quad_order);
        TMOPNewtonSolver solver(ir, solver_type);
        if (surface_fit_adapt) { solver.EnableAdaptiveSurfaceFitting(); }
        if (surface_fit_threshold > 0)
        {
          solver.SetTerminationWithMaxSurfaceFittingError(surface_fit_threshold);
        }
        // Provide all integration rules in case of a mixed mesh.
        solver.SetIntegrationRules(*irules, quad_order);
        if (solver_type == 0)
        {
          // Specify linear solver when we use a Newton-based solver.
          solver.SetPreconditioner(*S);
        }
        // For untangling, the solver will update the min det(T) values.
        solver.SetMinDetPtr(&min_detJ);
        solver.SetMaxIter(solver_iter);
        solver.SetRelTol(solver_rtol);
        solver.SetAbsTol(0.0);
        if (solver_art_type > 0)
        {
          solver.SetAdaptiveLinRtol(solver_art_type, 0.5, 0.9);
        }
        solver.SetPrintLevel(verbosity_level >= 1 ? 1 : -1);

        // hr-adaptivity solver.
        // If hr-adaptivity is disabled, r-adaptivity is done once using the
        // TMOPNewtonSolver.
        // Otherwise, "hr_iter" iterations of r-adaptivity are done followed by
        // "h_per_r_iter" iterations of h-adaptivity after each r-adaptivity.
        // The solver terminates if an h-adaptivity iteration does not modify
        // any element in the mesh.
        TMOPHRSolver hr_solver(*mesh, a, solver,
            x, move_bnd, hradaptivity,
            mesh_poly_deg, h_metric_id,
            n_hr_iter, n_h_iter);
        hr_solver.AddGridFunctionForUpdate(&x0);
        if (adapt_lim_const > 0.)
        {
          hr_solver.AddGridFunctionForUpdate(&adapt_lim_gf0);
          hr_solver.AddFESpaceForUpdate(&ind_fes);
        }
        hr_solver.Mult();


      }
    }

    if (o_mesh.is_curved() > 0) {
      auto cubic_curveVtk = oh::Mesh(o_mesh.comm()->library());
      cubic_curveVtk.set_comm(lib.world());
      build_cubic_curveVtk_3d(&o_mesh, &cubic_curveVtk, 5);
      std::string vtuPath = "/lore/joshia5/Meshes/curved/incl_aft_curveVtk.vtu";
      vtuPath += myid;
      oh::vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk, 2);
      auto cubic_wireframe = oh::Mesh(o_mesh.comm()->library());
      cubic_wireframe.set_comm(lib.world());
      build_cubic_wireframe_3d(&o_mesh, &cubic_wireframe, 5);
      vtuPath = "/lore/joshia5/Meshes/curved/incl_aft_wireframe.vtu";
      vtuPath += myid;
      oh::vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_wireframe, 1);
    }
    oh::vtk::write_parallel("after_adapt", &o_mesh);

    delete epsCoef;
    delete mesh;
    delete pmesh;
  } // end iterative adaptation loop

  return 0;
}

// The Permittivity is a required coefficient which may be defined in
// various ways so we'll determine the appropriate coefficient type here.
Coefficient *
SetupPermittivityCoefficient(int max_attr, // max attribute in the mesh
    double eps1, double eps2, // Permittivity of each phase
    const Array<int>& p1, const Array<int>& p2) // list of regions for each phase
{
   Vector epsilons(max_attr);
   epsilons = 0.0;
   for (int i=0; i<p1.Size(); i++)
     epsilons[p1[i]-1] = eps1;
   for (int i=0; i<p2.Size(); i++)
     epsilons[p2[i]-1] = eps2;
   Coefficient * coef = new PWConstCoefficient(epsilons);
   return coef;
}

// To produce a uniform electric field the potential can be set
// to (- Ex x - Ey y - Ez z).
double phi_bc_uniform(const Vector &x)
{
   double phi = 0.0;
   Vector euniform(3);
   euniform[0] = 0.;
   euniform[1] = 0.;
   euniform[2] = 1.;

   for (int i=0; i<x.Size(); i++)
   {
      phi -= x(i) * euniform(i);
   }

   return phi;
}
