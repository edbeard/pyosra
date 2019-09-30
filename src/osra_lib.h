/******************************************************************************
 OSRA: Optical Structure Recognition Application

 Created by Igor Filippov, 2007-2013 (igor.v.filippov@gmail.com)

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
 St, Fifth Floor, Boston, MA 02110-1301, USA
 *****************************************************************************/

// Header: osra_lib.h
//
// Defines types and functions of OSRA library.
//

#include <string> // std::string
#include <ostream> // std:ostream
#include <map>

//
// Section: Functions
//

// Function: osra_process_image()
//
// Parameters:
//      image_data - the binary image
//
// Returns:
//      0, if processing was completed successfully
int osra_process_image(
#ifdef OSRA_LIB
  const char *image_data,
  int image_length,
  std::ostream &structure_output_stream,
#else
  const std::string &input_file = "/home/edward/Documents/CSDE_Git/osra_install_pkg/osra-2.1.0-1/src/test.jpg",
  const std::string &output_file = "/home/edward/Documents/CSDE_Git/osra_install_pkg/osra-2.1.0-1/src/output",
#endif
  int rotate = 0,
  bool invert = false,
  int input_resolution = 0,
  double threshold = 0,
  int do_unpaper = 0,
  bool jaggy = false,
  bool adaptive_option = false,
  std::string output_format = "smi",
  std::string embedded_format = "",
  bool show_confidence = false,
  bool show_resolution_guess = false,
  bool show_page = false,
  bool show_coordinates = false,
  bool show_avg_bond_length = false,
  bool show_learning = false,
  const std::string &osra_dir = "/usr/local/bin",
  const std::string &spelling_file = "",
  const std::string &superatom_file = "",
  bool debug = false,
  bool verbose = false,
  const std::string &output_image_file_prefix = "",
  const std::string &resize = "",
  const std::string &preview = ""
);

// Instatiate and populate map and vector
std::vector<std::map<std::string, std::string> > initialize_rgroup();

/// Function for reading a single OSRA diagram
std::string read_diagram(
        const std::string &input_file = "/home/edward/cpp/hack_osra/src/test_rgroup.jpg",
        const char *image_data = "a",
        int image_length = 4,
        const std::string &output_file = "/home/edward/cpp/hack_osra/src/output",
        int rotate = 0,
        bool invert = false,
        int input_resolution = 0,
        double threshold = 0,
        int do_unpaper = 0,
        bool jaggy = false,
        bool adaptive_option = false,
        std::string output_format = "smi",
        std::string embedded_format = "",
        bool show_confidence = false,
        bool show_resolution_guess = false,
        bool show_page = false,
        bool show_coordinates = false,
        bool show_avg_bond_length = false,
        bool show_learning = false,
        const std::string &osra_dir = "/usr/local/bin",
        const std::string &spelling_file = "",
        const std::string &superatom_file = "",
        bool debug = false,
        bool verbose = true,
        const std::string &output_image_file_prefix = "",
        const std::string &resize = "",
        const std::string &preview = ""
);


std::vector<std::string> read_rgroup(
  std::vector<std::map<std::string, std::string> > list_of_rgroup_maps,
  const std::string &input_file = "/home/edward/cpp/hack_osra/src/test_rgroup.jpg",
  const char *image_data = "a",
  int image_length = 4,
  const std::string &output_file = "/home/edward/cpp/hack_osra/src/output",
  int rotate = 0,
  bool invert = false,
  int input_resolution = 0,
  double threshold = 0,
  int do_unpaper = 0,
  bool jaggy = false,
  bool adaptive_option = false,
  std::string output_format = "smi",
  std::string embedded_format = "",
  bool show_confidence = false,
  bool show_resolution_guess = false,
  bool show_page = false,
  bool show_coordinates = false,
  bool show_avg_bond_length = false,
  bool show_learning = false,
  const std::string &osra_dir = "/usr/local/bin",
  const std::string &spelling_file = "",
  const std::string &superatom_file = "",
  bool debug = false,
  bool verbose = false,
  const std::string &output_image_file_prefix = "",
  const std::string &resize = "",
  const std::string &preview = ""
);

void test_osra_lib(
#ifdef OSRA_LIB
  const std::string &output = "osra lib on",
#else
  const std::string &output = "osra lib off",
#endif
//  std::string output = "normal test",
  int pointless = 0
);