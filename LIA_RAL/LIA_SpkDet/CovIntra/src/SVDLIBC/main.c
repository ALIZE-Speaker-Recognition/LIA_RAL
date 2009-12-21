#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <strings.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "svdlib.h"

enum algorithms{LAS2};

/***********************************************************************
 *                                                                     *
 *                        timer()                                      *
 *            Returns elapsed cpu time (float, in seconds)             *
 *                                                                     *
 ***********************************************************************/
float timer(void) {
  long elapsed_time;
  struct rusage mytime;
  getrusage(RUSAGE_SELF,&mytime);
 
  /* convert elapsed time to milliseconds */
  elapsed_time = (mytime.ru_utime.tv_sec * 1000 +
                  mytime.ru_utime.tv_usec / 1000);
  
  /* return elapsed time in seconds */
  return((float)elapsed_time/1000.);
}

static long imin(long a, long b) {return (a < b) ? a : b;}

static void debug(char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
}

static void fatalError(char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  fprintf(stderr, "ERROR: ");
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\a\n");
  va_end(args);
  exit(1);
}

void printUsage(char *progname) {
  debug("SVD Version %s\n" 
        "written by Douglas Rohde based on code adapted from SVDPACKC\n\n", SVDVersion);
  debug("usage: %s [options] matrix_file\n", progname);
  debug("  -a algorithm   Sets the algorithm to use.  They include:\n"
        "       las2 (default)\n"
        "  -c infile outfile\n"
        "                 Convert a matrix file to a new format (using -r and -w)\n"
        "                 Then exit immediately\n"
        "  -d dimensions  Desired SVD triples (default is all)\n"
        "  -e bound       Minimum magnitude of wanted eigenvalues (1e-30)\n"
        "  -k kappa       Accuracy parameter for las2 (1e-6)\n"
        "  -i iterations  Algorithm iterations\n"
        "  -o file_root   Root of files in which to store resulting U,S,V\n"
        "  -r format      Input matrix file format\n"
        "       sth       SVDPACK Harwell-Boeing text format\n"
        "       st        Sparse text (default)\n"
        "       dt        Dense text\n"
        "       sb        Sparse binary\n"
        "       db        Dense binary\n"
        "  -v verbosity   Default 1.  0 for no feedback, 2 for more\n"
        "  -w format      Output matrix file format (see -r for formats)\n"
        "                   (default is dense text)\n");
  exit(1);
}


int main(int argc, char *argv[]) {
  extern char *optarg;
  extern int optind;
  int opt;

  SVDRec R = NULL;
  SMat A = NULL;

  char transpose = FALSE;
  int readFormat = SVD_F_ST;
  int writeFormat = SVD_F_DT;
  int algorithm = LAS2;
  int iterations = 0;
  int dimensions = 0;
  char *vectorFile = NULL;
  double las2end[2] = {-1.0e-30, 1.0e-30};
  double kappa = 1e-6;
  double exetime;

  while ((opt = getopt(argc, argv, "a:c:d:e:hk:i:o:r:tv:w:")) != -1) {
    switch (opt) {
    case 'a':
      if (!strcasecmp(optarg, "las2"))
        algorithm = LAS2;
      else fatalError("unknown algorithm: %s", optarg);
      break;
    case 'c':
      if (optind != argc - 1) printUsage(argv[0]);
      if (SVDVerbosity > 0) printf("Converting %s to %s\n", optarg, argv[optind]);
      if (SVD_IS_SPARSE(readFormat) && SVD_IS_SPARSE(writeFormat)) {
        SMat S = svdLoadSparseMatrix(optarg, readFormat);
        if (!S) fatalError("failed to read sparse matrix");
        if (transpose) {
          if (SVDVerbosity > 0) printf("  Transposing the matrix...\n");
          S = svdTransposeS(S); /* This leaks memory. */
        }
        svdWriteSparseMatrix(S, argv[optind], writeFormat);
      } else {
        DMat D = svdLoadDenseMatrix(optarg, readFormat);
        if (!D) fatalError("failed to read dense matrix");
        if (transpose) {
          if (SVDVerbosity > 0) printf("  Transposing the matrix...\n");
          D = svdTransposeD(D); /* This leaks memory. */
        }
        svdWriteDenseMatrix(D, argv[optind], writeFormat);
      }
      exit(0);
      break;
    case 'd':
      dimensions = atoi(optarg);
      if (dimensions < 0) fatalError("dimensions must be non-negative");
      break;
    case 'e':
      las2end[1] = atof(optarg);
      las2end[0] = -las2end[1];
      break;
    case 'h':
      printUsage(argv[0]);
      break;
    case 'k':
      kappa = atof(optarg);
      break;
    case 'i':
      iterations = atoi(optarg);
      break;
    case 'o':
      vectorFile = optarg;
      break;
    case 'r':
      if (!strcasecmp(optarg, "sth")) {
        readFormat = SVD_F_STH;
      } else if (!strcasecmp(optarg, "st")) {
        readFormat = SVD_F_ST;
      } else if (!strcasecmp(optarg, "dt")) {
        readFormat = SVD_F_DT;
      } else if (!strcasecmp(optarg, "sb")) {
        readFormat = SVD_F_SB;
      } else if (!strcasecmp(optarg, "db")) {
        readFormat = SVD_F_DB;
      } else fatalError("bad file format: %s", optarg);
      break;
    case 't':
      transpose = TRUE;
      break;
    case 'v':
      SVDVerbosity = atoi(optarg);
      /*if (SVDVerbosity) printf("Verbosity = %ld\n", SVDVerbosity);*/
      break;
    case 'w':
      if (!strcasecmp(optarg, "sth")) {
        writeFormat = SVD_F_STH;
      } else if (!strcasecmp(optarg, "st")) {
        writeFormat = SVD_F_ST;
      } else if (!strcasecmp(optarg, "dt")) {
        writeFormat = SVD_F_DT;
      } else if (!strcasecmp(optarg, "sb")) {
        writeFormat = SVD_F_SB;
      } else if (!strcasecmp(optarg, "db")) {
        writeFormat = SVD_F_DB;
      } else fatalError("bad file format: %s", optarg);
      break;
    default: printUsage(argv[0]);
    }
  }
  if (optind != argc - 1) printUsage(argv[0]);

  if (SVDVerbosity > 0) printf("Loading the matrix...\n");
  A = svdLoadSparseMatrix(argv[optind], readFormat);
  if (!A) fatalError("failed to read sparse matrix.  Did you specify the correct file type with the -r argument?");
  if (transpose) {
    if (SVDVerbosity > 0) printf("  Transposing the matrix...\n");
    SMat T = A;
    A = svdTransposeS(A);
    svdFreeSMat(T);
  }

  if (dimensions <= 0) dimensions = imin(A->rows, A->cols);

  exetime = timer();

  if (SVDVerbosity > 0) printf("Computing the SVD...\n");
  if (algorithm == LAS2) {
    if (!(R = svdLAS2(A, dimensions, iterations, las2end, kappa)))
      fatalError("error in svdLAS2");
  } else {
    fatalError("unknown algorithm");
  }

  exetime = timer() - exetime;
  if (SVDVerbosity > 0) {
    printf("\nELAPSED CPU TIME          = %6g sec.\n", exetime);
    printf("MULTIPLICATIONS BY A      = %6ld\n", 
           (SVDCount[SVD_MXV] - R->d) / 2 + R->d);
    printf("MULTIPLICATIONS BY A^T    = %6ld\n", 
           (SVDCount[SVD_MXV] - R->d) / 2);
  }

  if (vectorFile) {
    char filename[128];
    sprintf(filename, "%s-Ut", vectorFile);
    svdWriteDenseMatrix(R->Ut, filename, writeFormat);
    sprintf(filename, "%s-S", vectorFile);
    svdWriteDenseArray(R->S, R->d, filename, FALSE);
    sprintf(filename, "%s-Vt", vectorFile);
    svdWriteDenseMatrix(R->Vt, filename, writeFormat);
  }
  return 0;
}
