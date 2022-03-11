include(ExtraCmakeUtilities)
extra_find_package(METIS METIS METIS_DIR "include;Lib" "metis.h"
  "lib" "metis;metis4;metis5"
  "Paths to headers required by METIS." "Libraries required by METIS."
  CHECK_BUILD METIS_VERSION_5 FALSE
  "
#include <metis.h>
#include <cstddef> // So NULL is defined

int main()
{
    idx_t n = 10;
    idx_t nparts = 5;
    idx_t edgecut;
    idx_t* partitioning = new idx_t[10];
    idx_t* I = partitioning,
       * J = partitioning;

    idx_t ncon = 1;
    int err;
    idx_t options[40];

    METIS_SetDefaultOptions(options);
    options[10] = 1; // set METIS_OPTION_CONTIG

    err = METIS_PartGraphKway(&n,
                              &ncon,
                              I,
                              J,
                              (idx_t *) NULL,
                              (idx_t *) NULL,
                              (idx_t *) NULL,
                              &nparts,
                              (real_t *) NULL,
                              (real_t *) NULL,
                              options,
                              &edgecut,
                              partitioning);
    return err;
}
")
