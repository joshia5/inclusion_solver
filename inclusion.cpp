
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

  auto class_id = mesh->get_array<oh::LO> (0, "class_id");

  auto f = OMEGA_H_LAMBDA(oh::LO v) {
    auto h = oh::Vector<dim>();
    auto vtxError = error_c[v];
    for (oh::Int i = 0; i < dim; ++i) {
      h[i] = std::pow((error_des2/vtxError), 0.5)*length_c[v];
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
  opts.xfer_opts.type_map["zz_error"] = OMEGA_H_POINTWISE;
  opts.min_quality_allowed = 0.1;
  opts.max_length_allowed = 4.0*opts.max_length_desired;
  oh::Now t0 = oh::now();
  while (approach_metric(mesh, opts)) {
    printf("approach metric\n");
    adapt(mesh, opts);
    if (vtk_path) {
      writer.write();
    }
  }
  oh::Now t1 = oh::now();
  if (!myid) std::cout << "total time: " << (t1 - t0) << " seconds\n";

}

} // end anonymous namespace

// Permittivity Functions
Coefficient *
SetupPermittivityCoefficient(int max_attr, // max attribute in the mesh
    double eps1, double eps2, // Permittivity of each phase
    const Array<int>& p1, const Array<int>& p2); // list of regions for each phase

// Phi Boundary Condition
double phi_bc_uniform(const Vector &);

int main(int argc, char *argv[])
{
   MPI_Session mpi(argc, argv);
   int myid;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Read Omega_h mesh
     auto lib = oh::Library();
       oh::Mesh o_mesh(&lib);
  //
   oh::binary::read (
   "/lore/joshia5/Meshes/oh-mfem/inclusion_setup_5x5_2p.osh",
                    lib.world(), &o_mesh);

  double error_des = 0.0;

  //number of adaptation iterations
  int max_iter = 10;
  for (int Itr = 0; Itr < max_iter; Itr++)  {

    // problem constants and attribute and bdr attribute lists corresponding
    // to different regions and boundaries of the problem
    // NOTE: the list containing model tags are model dependent
    double kappa = 2.; // relative permittivity of phase 2 wrt vacuum
    //int num_substrate = 1; // number of regions in the substrate phase (1)
    //int num_inclusion = 1; // number of regions in the inclusion phase (2)
    //int substrate_regions[1] = {186};
    //int inclusion_regions[1] = {92};
    double epsilon1 = epsilon0_; // permittivity of substrate phase (1)
    double epsilon2 = kappa * epsilon0_; // permittivity of inclusion phase (2)
    bool isAmr = false;



    // Parse command-line options.
    const char *mesh_file = "";
    int order = 1;
    int maxit = 1;
    bool visualization = false;
    bool visit = false;

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
    args.AddOption(&maxit, "-maxit", "--max-amr-iterations",
        "Max number of iterations in the main AMR loop.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
        "--no-visualization",
        "Enable or disable GLVis visualization.");
    args.AddOption(&visit, "-visit", "--visit", "-no-visit", "--no-visit",
        "Enable or disable VisIt visualization.");
    args.Parse();
    if (!args.Good())
    {
      if (mpi.Root())
      {
        args.PrintUsage(cout);
      }
      return 1;
    }
    if (mpi.Root())
    {
      args.PrintOptions(cout);
    }

    // Read the (serial) mesh from the given mesh file on all processors.  We
    // can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
    // and volume meshes with the same code.
    //Mesh *mesh = new Mesh(mesh_file, 1, 1);
    ParMesh *pmesh = new ParOmegaMesh (MPI_COMM_WORLD, &o_mesh);

    int dim = pmesh->SpaceDimension();

    if (mpi.Root())
    {
      cout << "Starting initialization." << endl;
    }

    // Define a parallel mesh by a partitioning of the serial mesh. Refine
    // this mesh further in parallel to increase the resolution. Once the
    // parallel mesh is defined, the serial mesh can be deleted.

    pmesh->Finalize(true);


    // Create a coefficient describing the dielectric permittivity
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
    if (visualization)
    {
      Inclusion.InitializeGLVis();
    }

    // Initialize VisIt visualization
    VisItDataCollection visit_dc("Inclusion-AMR-Parallel", pmesh);

    if ( visit )
    {
      Inclusion.RegisterVisItFields(visit_dc);
    }
    if (mpi.Root()) { cout << "Initialization done." << endl; }

    // The main AMR loop. In each iteration we solve the problem on the current
    // mesh, visualize the solution, estimate the error on all elements, refine
    // the worst elements and update all objects to work with the new mesh.  We
    // refine until the maximum number of dofs in the nodal finite element space
    // reaches 10 million.
    const int max_dofs = 10000000;
    int it = 1;

      // Display the current number of DoFs in each finite element space
      Inclusion.PrintSizes();

      // Assemble all forms
      Inclusion.Assemble();

      // Solve the system and compute any auxiliary fields
      Inclusion.Solve();

      // Determine the current size of the linear system
      int prob_size = Inclusion.GetProblemSize();

      // Write fields to disk for VisIt
      if ( visit )
      {
        Inclusion.WriteVisItFields(it);
      }

      // Send the solution by socket to a GLVis server.
      if (visualization)
      {
        Inclusion.DisplayToGLVis();
      }


      // Check stopping criteria
      if (prob_size > max_dofs)
      {
        if (mpi.Root())
        {
          cout << "Reached maximum number of dofs, exiting..." << endl;
        }
        break;
      }

      // Wait for user input. Ask every 10th iteration.
      char c = 'c';
      if (mpi.Root() && (it % 10 == 0))
      {
        cout << "press (q)uit or (c)ontinue --> " << flush;
        cin >> c;
      }
      MPI_Bcast(&c, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

      if (c != 'c')
      {
        break;
      }

      Vector errors(pmesh->GetNE());
      Inclusion.GetErrorEstimates(errors);
      if (isAmr)
      {
        // Estimate element errors using the Zienkiewicz-Zhu error estimator.

        double local_max_err = errors.Max();
        double global_max_err;
        MPI_Allreduce(&local_max_err, &global_max_err, 1,
            MPI_DOUBLE, MPI_MAX, pmesh->GetComm());

        // Refine the elements whose error is larger than a fraction of the
        // maximum element error.
        const double frac = 0.1;
        double threshold = frac * global_max_err;
        if (mpi.Root()) { cout << "Refining ..." << endl; }
        pmesh->RefineByError(errors, threshold);
      }
      // adapt
      char Fname[128];
      sprintf(Fname,
          "inclusion_1x1_12k_4p.vtk");
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
      const double frac = 0.1;
      if (Itr == 0) error_des = frac*error_bar;

      printf("before adapt run case error max %1.10f error min %1.19f error_bar %1.10f error_des %1.10f\n",
          global_max_err, global_min_err, error_bar, error_des);
      //if (error_bar < error_des) break;

      stringstream ss;
      ss << "inclusion_iter_" << Itr;
      Inclusion.WriteToVtk(ss.str().c_str());

      if (global_max_err < error_des) {
        oh::vtk::write_parallel("before_break", &o_mesh);
        cout << "converged\n";
        break;
      }

      if (Itr < max_iter) run_case<3>(&o_mesh, Fname, Itr, myid, pOmesh, error_des);
      oh::vtk::write_parallel("after_adapt", &o_mesh);

      // Update the electrostatic solver to reflect the new state of the mesh.
      Inclusion.Update();

    delete epsCoef;
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
