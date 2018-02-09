#include "cpu_benchmark.hpp"
//#include "gpu_benchmark.hpp"

int main() {
    //subdiv_t<0,4,to_list_t<P<64, 64>, P<64, 1>>>::nothing;
    //subdiv_t<1, 4, to_list_t<P<64, 64>, P<64, 1>>>::nothing;
    //flip_t<1, subdiv_t<1, 4, to_list_t<P<64, 64>, P<64, 1>>>>::nothing;
    cpu::benchmark();
}