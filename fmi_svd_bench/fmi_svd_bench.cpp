#include <aws/core/Aws.h>
#include <aws/core/utils/json/JsonSerializer.h>
#include <aws/lambda-runtime/runtime.h>
#include <string>
#include <sstream>
#include <tcpunch.h>
#include <fcntl.h>
#include <Communicator.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <chrono>
#include "kiss.h"
#include "cblas.h"
#include "lapacke.h"
#include "mmio_dense.h"
#include "svd_utils.h"
#include "utils.h"
#include "fmi_wrapper.h"
#include "svd_algs_fmi.h"

using namespace aws::lambda_runtime;

static invocation_response my_handler(invocation_request const &req)
{
    using namespace Aws::Utils::Json;
    JsonValue json(req.payload);
    if (!json.WasParseSuccessful())
    {
        return invocation_response::failure("Failed to parse input JSON", "InvalidJSON");
    }

    auto v = json.View();
    std::string timestamp = v.GetString("timestamp");
    int nprocs = v.GetInteger("nprocs");
    int myrank = v.GetInteger("myrank");
    int m = v.GetInteger("rows");
    int n = v.GetInteger("cols");
    int p = v.GetInteger("ncomps");
    int seed = v.GetInteger("seed");

    double *Aloc, *Up, *Sp, *Vtp;

    auto comm = FMI::Communicator(myrank, nprocs, "fmi_svd_bench_config.json", timestamp);

    std::ostringstream res;

    auto t1 = std::chrono::high_resolution_clock::now();

    kiss_seed(myrank*seed);

    int s = n / nprocs;

    Aloc = (double*)malloc(m*s*sizeof(double));

    for (int i = 0; i < m*s; ++i)
    {
        Aloc[i] = 10*kiss_unirandf();
    }

    if (!myrank)
    {
        Up = (double*)malloc(m*p*sizeof(double));
        Sp = (double*)malloc(p*sizeof(double));
        Vtp = (double*)malloc(p*n*sizeof(double));
    }
    else
    {
        Up = Sp = Vtp = NULL;
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(t2-t1);
    res << myrank << " took " << elapsed.count() << " seconds to initialize\n";

    t1 = std::chrono::high_resolution_clock::now();
    if (svd_fmi(Aloc, Up, Sp, Vtp, m, n, p, 0, myrank, nprocs, comm) != 0)
    {
        return invocation_response::failure("svd_fmi failed", "svd_fmi_error");
    }
    t2 = std::chrono::high_resolution_clock::now();


    elapsed = std::chrono::duration_cast<std::chrono::seconds>(t2-t1);
    res << myrank << " took " << elapsed.count() << " seconds to run svd_fmi\n";

    if (!myrank)
    {
        free(Up);
        free(Sp);
        free(Vtp);
    }

    free(Aloc);

    return invocation_response::success(res.str(), "application/txt");
}


int main()
{
    using namespace Aws;

    auto handler_fn = [](aws::lambda_runtime::invocation_request const &req)
    {
        return my_handler(req);
    };

    run_handler(handler_fn);
    return 0;
}
