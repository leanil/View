////////////////////////////////////////////////////////////////////////////////
//
//  SYCL Sample Code
//
//  Author: Codeplay, April 2015
//
////////////////////////////////////////////////////////////////////////////////

/*!
Please note that this project will be automatically generated a Win32 project,
if you are using a 64bit ComputeCpp installation please switch to the x64
configuration before building.
*/

/*!
If you encounter an application cannot start because SYCL.dll or SYCL_d.dll is
missing error, then you may have to restart your computer in order to update the
environment variable COMPUTECPPROOT.
*/

#include "SYCL/sycl.hpp"
#include <vector>
#include <iostream>

using namespace cl::sycl;

static const int n = 1 << 20;

int gpu_benchmark() {

    /* Create 1024 element arrays for input data and output. */
    std::vector<std::vector<float>> data(3);
    for (auto& v : data) { v.resize(n); }

    /* Initialize input data with values and output data with zeroes. */
    for (int i = 0; i < n; i++) {
        data[0][i] = (float)i;
        data[1][i] = (float)(n - i);
        data[2][i] = 0.0f;
    }

    /* Wrap all SYCL structures and function calls with a try-catch block to catch
    * SYCL exceptions. */
    try {

        /* Create a default_selector to select a device to execute on. */
        cpu_selector mySelector;

        /* Create a queue from the default_selector to create an implicit context
        * and queue to execute with. */
        property_list plist{ cl::sycl::property::queue::enable_profiling {} };
        queue myQueue(mySelector, plist);

        /* Create a scope to control data synchronisation of buffer objects. */
        {
            /* Create 1 dimensionsal buffers for the input and output data of size
            * 1024. */
            buffer<float, 1> inputBufferA(data[0].data(), range<1>(n));
            buffer<float, 1> inputBufferB(data[1].data(), range<1>(n));
            buffer<float, 1> outputBuffer(data[2].data(), range<1>(n));

            /* Submit a command_group to execute from the queue. */
            auto e = myQueue.submit([&](handler &cgh) {

                /* Create accessors for accessing the input and output data within the
                * kernel. */
                auto inputPtrA = inputBufferA.get_access<access::mode::read>(cgh);
                auto inputPtrB = inputBufferB.get_access<access::mode::read>(cgh);
                auto outputPtr = outputBuffer.get_access<access::mode::write>(cgh);

                /* Enqueue a kernel called 'vector_add', with a global work size of {
                * 16, 8, 8 } and a local work size of { 4, 2, 2 }. */
                cgh.parallel_for<class vector_add>(range<1>(n), [=](item<1> item) {

                    /* Retreive the linear global id for the current work item. */
                    size_t idx = item.get_linear_id();

                    /* Use the linear global id to add the respective element of each
                    * input accessor together and assign them to the respective
                    * element of the output accessor. */
                    outputPtr[idx] = inputPtrA[idx] + inputPtrB[idx];
                });
            });
            myQueue.wait();
            cl::sycl::cl_ulong t0, t1;
            auto status1 = clGetEventProfilingInfo(e.get(), CL_PROFILING_COMMAND_START, sizeof(cl::sycl::cl_ulong), &t0, nullptr);
            auto status2 = clGetEventProfilingInfo(e.get(), CL_PROFILING_COMMAND_END, sizeof(cl::sycl::cl_ulong), &t1, nullptr);
            if (status1 == CL_SUCCESS && status2 == CL_SUCCESS)
            {
                printf("time: %f ms\n", (t1 - t0) / 1'000'000.0f);
            }
            else
            {
                std::cout << "Could not get profiling info\n";
            }
        }

    }
    catch (exception e) {

        /* In the case of an exception being throw, print theerror message and
        * return 1. */
        std::cout << e.what();
        return 1;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////