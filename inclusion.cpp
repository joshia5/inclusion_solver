
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

namespace oh = Omega_h;

namespace { // anonymous namespace

/* parts of this function is derived from the file
 * ugawg_linear.cpp of omega_h source code
 */
template <oh::Int dim>
static void set_target_metric(oh::Mesh* mesh, oh::Int scale, ParOmegaMesh
  *pOmesh) {
  auto coords = mesh->coords();
  auto target_metrics_w = oh::Write<oh::Real>
    (mesh->nverts() * oh::symm_ncomps(dim));
  pOmesh->ProjectFieldElementtoVertex (mesh, "zz_error");
  auto zz_error = mesh->get_array<oh::Real> (0, "zz_error");
  auto f = OMEGA_H_LAMBDA(oh::LO v) {
    auto x = coords[v*dim];
    auto y = coords[v*dim + 1];
    auto z = coords[v*dim + 2];
    auto h = oh::Vector<dim>();
    auto vtxError = zz_error[v];
    for (oh::Int i = 0; i < dim; ++i)
      h[i] = 0.004/std::pow(std::abs(vtxError), 0.6);
      //h[i] = 0.001/std::pow(std::abs(vtxError), 0.6);//1k, 0.33mil
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  oh::parallel_for(mesh->nverts(), f);
  mesh->set_tag(oh::VERT, "target_metric", oh::Reals(target_metrics_w));
}

/* parts of this function is derived from the file
 * ugawg_linear.cpp of omega_h source code
 */
template <oh::Int dim>
void run_case(oh::Mesh* mesh, char const* vtk_path, oh::Int scale,
              const oh::Int myid, ParOmegaMesh *pOmesh) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = get_implied_metrics(mesh);
  mesh->add_tag(oh::VERT, "metric", oh::symm_ncomps(dim), implied_metrics);
  mesh->add_tag<oh::Real>(oh::VERT, "target_metric", oh::symm_ncomps(dim));
  set_target_metric<dim>(mesh, scale, pOmesh);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->ask_lengths();
  mesh->ask_qualities();
  oh::vtk::FullWriter writer;
  if (vtk_path) {
    writer = oh::vtk::FullWriter(vtk_path, mesh);
    writer.write();
  }
  auto opts = oh::AdaptOpts(mesh);
  opts.verbosity = oh::EXTRA_STATS;
  opts.length_histogram_max = 2.0;
  opts.max_length_allowed = opts.max_length_desired * 4.0;
  opts.min_quality_allowed = 0.00001;
  opts.xfer_opts.type_map["zz_error"] = OMEGA_H_POINTWISE;
  oh::Now t0 = oh::now();
  while (approach_metric(mesh, opts)) {
    adapt(mesh, opts);
    if (mesh->has_tag(oh::VERT, "target_metric")) set_target_metric<dim>(mesh,
                      scale, pOmesh);
    if (vtk_path) writer.write();
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
   oh::binary::read ("/lore/joshia5/Meshes/oh-mfem/inclusion_1x1_12k_4p.osh",
                    lib.world(), &o_mesh);

  //number of adaptation iterations
  int max_iter = 1;
  for (int Itr = 0; Itr < max_iter; Itr++)  {

   // problem constants and attribute and bdr attribute lists corresponding
   // to different regions and boundaries of the problem
   // NOTE: the list containing model tags are model dependent
   double kappa = 2.; // relative permittivity of phase 2 wrt vacuum
   int num_substrate = 1; // number of regions in the substrate phase (1)
   int num_inclusion = 1; // number of regions in the inclusion phase (2)
   int substrate_regions[1] = {186};
   int inclusion_regions[1] = {92};
   double epsilon1 = epsilon0_; // permittivity of substrate phase (1)
   double epsilon2 = kappa * epsilon0_; // permittivity of inclusion phase (2)
   bool isAmr = false;



   // Parse command-line options.
   const char *mesh_file = "";
   int order = 1;
   int maxit = 20;
   bool visualization = true;
   bool visit = true;

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

   int sdim = pmesh->SpaceDimension();

   if (mpi.Root())
   {
      cout << "Starting initialization." << endl;
   }


   // Define a parallel mesh by a partitioning of the serial mesh. Refine
   // this mesh further in parallel to increase the resolution. Once the
   // parallel mesh is defined, the serial mesh can be deleted.
   //ParMesh pmesh(MPI_COMM_WORLD, *mesh);
   //delete mesh;

   // Make sure tet-only meshes are marked for local refinement.
   //pmesh.Finalize(true);


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
   for (int it = 1; it <= maxit; it++)
   {
      if (mpi.Root())
      {
         cout << "\nAMR Iteration " << it << endl;
      }

      // Display the current number of DoFs in each finite element space
      Inclusion.PrintSizes();

      // Assemble all forms
      Inclusion.Assemble();

      // Solve the system and compute any auxiliary fields
      Inclusion.Solve();

      stringstream ss;
      ss << "inclusion_iter_" << it;
      Inclusion.WriteToVtk(ss.str().c_str());
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

      if (mpi.Root())
      {
         cout << "AMR iteration " << it << " complete." << endl;
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
      if (it == maxit)
      {
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

      if (isAmr)
      {
	// Estimate element errors using the Zienkiewicz-Zhu error estimator.
	Vector errors(pmesh->GetNE());
	Inclusion.GetErrorEstimates(errors);

	double local_max_err = errors.Max();
	double global_max_err;
	MPI_Allreduce(&local_max_err, &global_max_err, 1,
		      MPI_DOUBLE, MPI_MAX, pmesh->GetComm());

	// Refine the elements whose error is larger than a fraction of the
	// maximum element error.
	const double frac = 0.7;
	double threshold = frac * global_max_err;
	if (mpi.Root()) { cout << "Refining ..." << endl; }
	pmesh->RefineByError(errors, threshold);
      }

      // Update the electrostatic solver to reflect the new state of the mesh.
      Inclusion.Update();

      if (pmesh->Nonconforming() && mpi.WorldSize() > 1)
      {
         if (mpi.Root()) { cout << "Rebalancing ..." << endl; }
         pmesh->Rebalance();

         // Update again after rebalancing
         Inclusion.Update();
      }
   }

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
