#include "nekAscent.hpp"

MPI_Comm comm;
ascent::Ascent mAscent;

static double *xm1, *ym1, *zm1; // ptr to mesh

// work arrays
static double *x, *y, *z; // mesh
static hlong *con;        // connectivity
static double *flds;      // solution fields

static int nfld;          // # solutions
static hlong ntot, ltot;
static fields userFieldList;
static conduit::Node mesh_data;

static bool if3d = false;
static bool initCalled = false;
static bool setupCalled = false;
static bool updateMesh = true;

static int rank;
static int verbose = 1;

constexpr unsigned int NSCALAR_MAX = 100;

static void copy(double *a, const double *b, const hlong n) {
  for (hlong i=0;i<n;i++) a[i] = b[i];
}
static void rzero(double *a, const hlong n) {
  for (hlong i=0;i<n;i++) a[i] = 0.0;
}

static const std::string scalarDigitStr(int i)
{ 
  const int scalarWidth = std::to_string(NSCALAR_MAX - 1).length();
  std::stringstream ss;
  ss << std::setfill('0') << std::setw(scalarWidth) << i;
  return ss.str();
};


static void ascent_init() {

  const double tStart = MPI_Wtime();
  if (rank == 0) {
    printf("Initialize Ascent ...\n");
    fflush(stdout);
  }

  conduit::Node ascent_opts;
  ascent_opts["mpi_comm"]=MPI_Comm_c2f(comm);

//  ascent_opts["runtime/type"] = "ascent";
/*
  std::string backend; // TODO: do nothing? openmp? dpcpp? support mismatch?
  platform->options.getArgs("THREAD MODEL", backend);
  if (backend == "CUDA") {
    ascent_opts["runtime/backend"] = "cuda";
  } else if (backend == "HIP") {
    ascent_opts["runtime/backend"] = "hip";
  } else if (backend == "CPU" || backend == "SERIAL") {
    ascent_opts["runtime/backend"] = "serial";
  } else {
    ascent_opts["runtime/backend"] = "openmp";
  }
*/
  if (verbose) { // FIXME: This doesn't do much?
    ascent_opts["ascent_info"] = "verbose";
    ascent_opts["messages"] = "verbose";
  }
  mAscent.open(ascent_opts);

  const double tSetup = MPI_Wtime() - tStart; 
  if (rank == 0) {
    printf("done (%gs)\n\n", tSetup);
    std::cout << ascent::about() << std::endl; 
    fflush(stdout);
  }

  initCalled = true;
}


void nekascent_setup(MPI_Comm comm_in,
                     int ndim_, int lx1_, int ly1_, int lz1_, 
                     int nelt_, int lelt_, int nfldt_,
                     double *xm1_, double *ym1_, double *zm1_,
                     double *vx, double *vy, double *vz, 
                     double *pm1, double *t) {

  // mpi
  comm = comm_in;
  MPI_Comm_rank(comm, &rank);

  // verbose level
  if (rank == 0) {
    const char *val = getenv("NEKASCENT_VERBOSE_LEVEL");
    if (val != NULL) verbose = atoi(val);
  }
  MPI_Bcast(&verbose, 1, MPI_INT, 0, comm);

  // init
  if (!initCalled) ascent_init();

  // timer init
  const double tStart = MPI_Wtime();
  if (rank == 0) {
    printf("Setup Ascent fields...");
    fflush(stdout);
  }

  const int ndim = ndim_, lx1 = lx1_, ly1 = ly1_, lz1 = lz1_;
  const int nelt = nelt_, lelt = lelt_, nfldt = nfldt_;

  xm1 = xm1_;
  ym1 = ym1_;
  zm1 = zm1_;

  if3d = (ndim==3) ? true : false;
  nfld = ndim + 1 + nfldt; // vel + pr + scalars
  ntot = lx1*ly1*lz1*nelt;
  ltot = lx1*ly1*lz1*lelt;

  hlong nvtx, ncells;
  if (if3d) {
    ncells = (lx1-1)*(ly1-1)*(lz1-1)*nelt;
    nvtx = ncells*8;
  } else {
    ncells = (lx1-1)*(ly1-1)*nelt;
    nvtx = ncells*4;
  }

  if (rank==0)// TODO: verbose
    printf("setup: %lld %lld %lld %d\n",ntot,ltot,nvtx,nfld);

  // allocate work arrays
  x = (double *) malloc(sizeof(double)*ltot);
  y = (double *) malloc(sizeof(double)*ltot);
  if (if3d) z = (double *) malloc(sizeof(double)*ltot);
  flds = (double *) malloc(sizeof(double)*ltot*nfld);
  con = (hlong *) malloc(sizeof(hlong)*nvtx);

  // put data pointer to fields arrays
  userFieldList.push_back(std::make_tuple("vel_x", vx));
  userFieldList.push_back(std::make_tuple("vel_y", vy));
  if (if3d) userFieldList.push_back(std::make_tuple("vel_z", vz));
  userFieldList.push_back(std::make_tuple("pressure", pm1));
  for (int is=0; is<nfldt; is++) {
    const std::string sid = (is==0) ? "temperature" : "scalar" + scalarDigitStr(is);
    double *tmpFld = t + is*ltot;
    userFieldList.push_back(std::make_tuple(sid, tmpFld));
  }

  // calculate connectivity
  if (if3d) {
    hlong id = 0;
    for(int ie=0; ie<nelt; ie++)
      for(int iz=0; iz<lz1-1; iz++)
        for(int iy=0; iy<ly1-1; iy++)
          for(int ix=0; ix<lx1-1; ix++) {
            con[id+0] = ((ie * lz1 + iz) * ly1 + iy) * lx1 + ix;
            con[id+1] = con[id+0] + 1;
            con[id+2] = con[id+0] + lx1+1;
            con[id+3] = con[id+0] + lx1;
            con[id+4] = con[id+0] + lx1*ly1;
            con[id+5] = con[id+1] + lx1*ly1;
            con[id+6] = con[id+2] + lx1*ly1;
            con[id+7] = con[id+3] + lx1*ly1;
            id += 8;  
          }
  } else {
    hlong id = 0;
    for(int ie=0; ie<nelt; ie++)
      for(int iy=0; iy<ly1-1; iy++)
        for(int ix=0; ix<lx1-1; ix++) {
          con[id+0] = (ie * ly1 + iy) * lx1 + ix;
          con[id+1] = con[id+0] + 1;
          con[id+2] = con[id+0] + lx1+1;
          con[id+3] = con[id+0] + lx1;
          id += 4;
        }
  }

  // setup ascent data pointer
  mesh_data["coordsets/coords/type"] = "explicit";
  mesh_data["coordsets/coords/values/x"].set_external(x, ntot);
  mesh_data["coordsets/coords/values/y"].set_external(y, ntot);
  if (if3d) mesh_data["coordsets/coords/values/z"].set_external(z, ntot);

  mesh_data["topologies/mesh/type"]           = "unstructured";
  mesh_data["topologies/mesh/coordset"]       = "coords";
  if (if3d) {
    mesh_data["topologies/mesh/elements/shape"] = "hex"; 
  } else {
    mesh_data["topologies/mesh/elements/shape"] = "quad"; // TODO: need test
  }
  mesh_data["topologies/mesh/elements/connectivity"].set_external(con, nvtx);

  // fields
  int ifld = 0;
  for(auto& entry : userFieldList) {
    std::string fieldName = std::get<0>(entry);
    double *tmpFld = flds + ifld * ltot;
    mesh_data["fields/" + fieldName + "/association"]  = "vertex";
    mesh_data["fields/" + fieldName + "/topology"]     = "mesh";
    mesh_data["fields/" + fieldName + "/values"].set_external(tmpFld, ntot);
    ifld++;
  }

  const double tSetup = MPI_Wtime() - tStart;  // TODO: add Nfields and size?
  if (rank == 0) {
    printf("done (%gs)\n\n", tSetup);
    fflush(stdout);
  }

  setupCalled = true;
}


void nekascent_update(double time, int istep, int movingMesh) {

  if (!setupCalled) {
    printf("ascent_update called prior to ascent_setup!\n");
    exit(EXIT_FAILURE);
  }

  const double tStart = MPI_Wtime();

  // copy data
  mesh_data["state/cycle"] = istep;
  mesh_data["state/time"] = time;

  if (updateMesh || movingMesh) { // update mesh at first call or moving mesh
    copy(x, xm1, ntot);
    copy(y, ym1, ntot);
    if (if3d) copy(z, zm1, ntot);
    updateMesh = false;
    // todo add verbose print
  }

  rzero(flds, nfld*ltot);
  int ifld = 0;
  for(auto& entry : userFieldList) {
    double *dataPtr = std::get<1>(entry);
    double *tmpFld = flds + ifld*ltot;
    copy(tmpFld, dataPtr, ntot);
    ifld++;
  }

  conduit::Node verify_info;
  if (!conduit::blueprint::mesh::verify(mesh_data, verify_info)) {
    if (rank == 0)
      CONDUIT_INFO("blueprint verify failed!" + verify_info.to_yaml());
  } else {
/*
  // TODO:
  std::string actions_yaml = R"ST(
-       
  action: "add_queries"
  queries:
    q1:
      params:
        expression: "cycle()"
        name: "cycle"
-
  action: "add_triggers"
  triggers:
    t1:
      params:
        condition: "cycle % 10 == 0"
        actions_file : "trigger_render.yaml"
)ST";
*/

    mAscent.publish(mesh_data);
    conduit::Node actions;
//    actions.parse(actions_yaml,"yaml");
    mAscent.execute(actions);
  }

  const double tTotal = MPI_Wtime() - tStart;  
  if (rank == 0) {
    printf("Ascent update at step=%d time=%g, (%gs)\n", istep, time, tTotal);
    fflush(stdout);
  }

}


void nekascent_finalize() {
  mAscent.close();
  free(x);
  free(y);
  free(z);
  free(flds);
  free(con);
}



void fnekascent_setup(int *comm_in,
                      int *ndim_, int *lx1_, int *ly1_, int *lz1_,
                      int *nelt_, int *lelt_, int *nfldt_,
                      double *xm1_, double *ym1_, double *zm1_,
                      double *vx, double *vy, double *vz,
                      double *pm1, double *t) {

  MPI_Comm comm_ = MPI_Comm_f2c(*comm_in);
  nekascent_setup(comm_,
                  *ndim_, *lx1_, *ly1_, *lz1_,
                  *nelt_, *lelt_, *nfldt_,
                  xm1_, ym1_, zm1_,
                  vx, vy, vz,
                  pm1, t);
}

void fnekascent_update(double *time, int *istep, int *movingMesh) {
  nekascent_update(*time, *istep, *movingMesh);
}

void fnekascent_finalize() {
  nekascent_finalize();
}
