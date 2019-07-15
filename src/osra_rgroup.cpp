#include "iostream"
#include "stdio.h"
#include "osra_common.h"
#include "osra_lib.h"
#include "string"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;

//int main()
//{
////    std::cout << "file run" << std::endl;
////    std::cout << distance(3, 3,3,4) << std::endl;
////
////
////    const std::string input = "/home/edward/cpp/hack_osra/src/test_rgroup.jpg";
////    const std::string output = "/home/edward/cpp/hack_osra/src/output";
////
////    std::vector<std::string> var_vec{"R"};
////    std::vector<std::string> val_vec{"CH3"};
////    std::vector<std::map<std::string, std::string> > dummy_input;
//
//   // int smiles = hack_osra_process_image(dummy_input);
//    //std::cout << smiles << std::endl;
//
//    return 0;
//}

PYBIND11_MODULE(osra_rgroup, m){

    m.doc() = "pybind test example";
    // Function for resolving R-groups
    m.def("read_rgroup", &read_rgroup, "Extension of OSRA to resolve generic chemical diagrams (employed by Chemical Schematic Diagram Extractor)",
    py::arg("list_of_rgroup_maps"),
    py::arg("input_file"),
    py::arg("image_data") = "a",
    py::arg("image_length") = 4,
    py::arg("output_file") = "/home/edward/cpp/hack_osra/src/output",
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
    py::arg("verbose") = false,
    py::arg("output_image_file_prefix") = "",
    py::arg("resize") = "",
    py::arg("preview") = ""
    );

    // Function for directly resolving images
    m.def("read_diagram", &read_diagram, "Python wrapper of OSRA for resolving chemical diagrams",
    py::arg("input_file"),
    py::arg("image_data") = "a",
    py::arg("image_length") = 4,
    py::arg("output_file") = "/home/edward/cpp/hack_osra/src/output",
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
    py::arg("verbose") = false,
    py::arg("output_image_file_prefix") = "",
    py::arg("resize") = "",
    py::arg("preview") = ""
);
}