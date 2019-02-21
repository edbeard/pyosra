#include "iostream"
#include "stdio.h"
#include "osra_common.h"
#include "osra_lib.h"
#include "string"
#include <pybind11/pybind11.h>

namespace py = pybind11;

int main()
{
    std::cout << "file run" << std::endl;
    std::cout << distance(3, 3,3,4) << std::endl;

   
    const std::string input = "/home/edward/Documents/CSDE_Git/osra_install_pkg/osra-2.1.0-1/src/test.jpg";
    const std::string output = "/home/edward/Documents/CSDE_Git/osra_install_pkg/osra-2.1.0-1/src/output";

    int smiles = hack_osra_process_image("a", 4);
    std::cout << smiles << std::endl;

    return 0;
}

PYBIND11_MODULE(test, m){
    m.doc() = "pybind test example";
    m.def("hack_osra_process_image", &hack_osra_process_image, "Hack of OSRA to create smiles for input and output",
    py::arg("image_data") = "a",
    py::arg("image_length") = 4,
    py::arg("input_file") = "/home/edward/Documents/CSDE_Git/osra_install_pkg/osra-2.1.0-1/src/test.jpg",
    py::arg("output_file") = "/home/edward/Documents/CSDE_Git/osra_install_pkg/osra-2.1.0-1/src/output",
    py::arg("rotate") = 0,
    py::arg("invert") = false,
    py::arg("input_resolution") = 0,
    py::arg("threshold") = 0.,
    py::arg("do_unpaper") = 0,
    py::arg("jaggy") = false,
    py::arg("adaptive_option") = false,
    py::arg("output_format") = "smi",
    py::arg("embedded_format") = "",
    py::arg("show_confidence") = false,
    py::arg("show_resolution_guess") = false,
    py::arg("show_page") = false,
    py::arg("show_coordinates") = false,
    py::arg("show_avg_bond_length") = false,
    py::arg("show_learning") = false,
    py::arg("osra_dir") = "/usr/local/bin",
    py::arg("spelling_file") = "",
    py::arg("superatom_file") = "",
    py::arg("debug") = false,
    py::arg("verbose") = true,
    py::arg("output_image_file_prefix") = "",
    py::arg("resize") = "",
    py::arg("preview") = ""
    );

}