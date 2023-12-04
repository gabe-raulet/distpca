#include <aws/lambda-runtime/runtime.h>

using namespace aws::lambda_runtime;

invocation_response my_handler(invocation_request const& req)
{
    if (req.payload.length() > 42)
    {
        return invocation_response::failure("FAILURE", "BAD TYPE OF FAILURE");
    }

    return invocation_response::success("json payload here", "application/json");
}

int main()
{
   run_handler(my_handler);
   return 0;
}
