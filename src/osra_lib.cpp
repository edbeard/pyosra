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

#include <stdio.h> // fclose
#include <stdlib.h> // malloc(), free()
#include <math.h> // fabs(double)
#include <float.h> // FLT_MAX
#include <limits.h> // INT_MAX

#include <list> // sdt::list
#include <vector> // std::vector
#include <set>
#include <algorithm> // std::sort, std::min(double, double), std::max(double, double)
#include <iostream> // std::ostream, std::cout
#include <fstream> // std::ofstream, std::ifstream
#include <sstream> // std:ostringstream

#include <Magick++.h>

extern "C" {
#include <potracelib.h>
#include <pgm2asc.h>
}

#include <openbabel/oberror.h>
#include <poppler/cpp/poppler-document.h>
#include <poppler/cpp/poppler-page-renderer.h>

#include "osra.h"
#include "osra_grayscale.h"
#include "osra_segment.h"
#include "osra_fragments.h"
#include "osra_labels.h"
#include "osra_thin.h"
#include "osra_common.h"
#include "osra_structure.h"
#include "osra_lib.h"
#include "osra_ocr.h"
#include "osra_openbabel.h"
#include "osra_reaction.h"
#include "osra_anisotropic.h"
#include "osra_stl.h"
#include "unpaper.h"
#include "config.h" // DATA_DIR

using namespace Magick;

void set_select_resolution(std::vector<int>  &select_resolution, int input_resolution)
{
  if (input_resolution == 0)
    {
      select_resolution[0] = 72;
      select_resolution[1] = 150;
      select_resolution[2] = 300;
      select_resolution[3] = 300;  // No thinning
      select_resolution[4] = 500;
    }
}

double set_threshold(double threshold,int resolution)
{
  double THRESHOLD_BOND = threshold;

  if (THRESHOLD_BOND < 0.0001)
    {
      if (resolution >= 150)
        {
          THRESHOLD_BOND = THRESHOLD_GLOBAL;
        }
      else
        {
          THRESHOLD_BOND = THRESHOLD_LOW_RES;
        }
    }
  return THRESHOLD_BOND;
}

int load_superatom_spelling_maps(
    std::map<std::string, std::string> &spelling, std::map<std::string, std::string> &superatom,
    const std::string &osra_dir, const std::string &spelling_file,
    const std::string &superatom_file, bool verbose)
{
// Loading the program data files into maps:
  if (!((spelling_file.length() != 0 && load_config_map(spelling_file, spelling))
        || load_config_map(std::string(DATA_DIR) + "/" + SPELLING_TXT, spelling) || load_config_map(osra_dir + "/" + SPELLING_TXT, spelling)))
    {
      std::cerr << "Cannot open " << SPELLING_TXT << " file (tried locations \"" << DATA_DIR << "\", \"" << osra_dir
                << "\"). Specify the custom file location via -l option." << std::endl;
      return ERROR_SPELLING_FILE_IS_MISSING;
    }

  if (!((superatom_file.length() != 0 && load_config_map(superatom_file, superatom))
        || load_config_map(std::string(DATA_DIR) + "/" + SUPERATOM_TXT, superatom) || load_config_map(osra_dir + "/"
            + SUPERATOM_TXT, superatom)))
    {
      std::cerr << "Cannot open " << SUPERATOM_TXT << " file (tried locations \"" << DATA_DIR << "\", \"" << osra_dir
                << "\"). Specify the custom file location via -a option." << std::endl;
      return ERROR_SUPERATOM_FILE_IS_MISSING;
    }

  if (verbose)
    std::cout << "spelling (size: " << spelling.size() << ") and superatom (size: " << superatom.size() << ") dictionaries are loaded." << std::endl;
  return 0;
}

void create_thick_box(Image &orig_box,Image &thick_box,int &width,int &height,int &resolution,int &working_resolution,double &box_scale,
                      ColorGray bgColor, double THRESHOLD_BOND, int res_iter, bool &thick, bool jaggy)
{
  if (resolution >= 300)
    {
      int max_hist;
      double nf45;
      double nf =
        noise_factor(orig_box, width, height, bgColor, THRESHOLD_BOND, resolution, max_hist, nf45);

      //if (max_hist < 5) thick = false;
      if (res_iter == NUM_RESOLUTIONS-2)  // no thinning
	{
	  working_resolution = 300;
	  thick_box = orig_box;
	  width = thick_box.columns();
	  height = thick_box.rows();
	  thick = false;
	  //nf = noise_factor(orig_box, width, height, bgColor, THRESHOLD_BOND, resolution, max_hist, nf45);
	}
      if (res_iter == NUM_RESOLUTIONS-1)
        {
          if (max_hist >= 6)
            {
              int new_resolution = max_hist * 300 / 4;
              int percent = (100 * 300) / new_resolution;
              resolution = new_resolution;
              std::ostringstream scale;
              scale << percent << "%";
              orig_box.scale(scale.str());
              box_scale /= (double) percent/100;
              working_resolution = 300;
              thick_box = orig_box;
              width = thick_box.columns();
              height = thick_box.rows();
	      nf = noise_factor(orig_box, width, height, bgColor, THRESHOLD_BOND, resolution, max_hist, nf45);
            }
          else
            {
              resolution = 500;
              int percent = (100 * 300) / resolution;
              std::ostringstream scale;
              scale << percent << "%";
              orig_box.scale(scale.str());
              box_scale /= (double) percent/100;
              working_resolution = 300;
              thick_box = orig_box;
              width = thick_box.columns();
              height = thick_box.rows();
              thick = false;
              nf = noise_factor(orig_box, width, height, bgColor, THRESHOLD_BOND, resolution,  max_hist, nf45);
            }
        }
      if (jaggy)
        {
          orig_box.scale("50%");
          box_scale *= 2;
          working_resolution = 150;
          //orig_box.scale("33%");
          //box_scale *= 3;
          // working_resolution = 100;
          thick_box = orig_box;
          width = thick_box.columns();
          height = thick_box.rows();
        }
      else if (nf > 0.5 && nf < 1. && max_hist <= 6)// && res_iter != 3 && max_hist <= 6)
        try
          {
            thick_box = anisotropic_smoothing(orig_box, width, height, 20, 0.3, 1.0, 0.6, 2);
          }
        catch (...)
          {
            thick_box = orig_box;
          }
      /*else if (nf45 > 0.9 && nf45 < 1.2 && max_hist == 3)
      {
      //orig_box = anisotropic_smoothing(thick_box, width, height, 60, 0.3, 0.6, 4., 2.);
      orig_box.scale("50%");
      thick_box = orig_box;
      //working_resolution = 150;
      width = thick_box.columns();
      height = thick_box.rows();
      //thick = false;
      }*/
      else
        thick_box = orig_box;
    }
  else if (resolution < 300 && resolution > 150)
    {
      int nw = width * 300 / resolution;
      int nh = height * 300 / resolution;
      thick_box = anisotropic_scaling(orig_box, width, height, nw, nh);
      width = thick_box.columns();
      height = thick_box.rows();
      int percent = (100 * 300) / resolution;
      std::ostringstream scale;
      scale << percent << "%";
      orig_box.scale(scale.str());
      box_scale /= (double) percent/100;
      working_resolution = 300;
    }
  else
    thick_box = orig_box;
}

potrace_state_t * const  raster_to_vector(Image &box,ColorGray bgColor, double THRESHOLD_BOND,int width,int height,int working_resolution)
{
  potrace_param_t * const param = potrace_param_default();
  param->alphamax = 5e-324; // this has been changed in potrace-1.11 
  //param->turnpolicy = POTRACE_TURNPOLICY_MINORITY;
  param->turdsize = 0;

  param->turnpolicy = POTRACE_TURNPOLICY_MINORITY;
  double c_width = 1. * width * 72 / working_resolution;
  double c_height = 1. * height * 72 / working_resolution;
  if (c_height * c_width < SMALL_PICTURE_AREA)
    param->turnpolicy = POTRACE_TURNPOLICY_BLACK;

  potrace_bitmap_t * const bm = bm_new(width, height);
  for (int i = 0; i < width; i++)
    for (int j = 0; j < height; j++)
      BM_PUT(bm, i, j, get_pixel(box, bgColor, i, j, THRESHOLD_BOND));

  potrace_state_t * const st = potrace_trace(param, bm);
  if (bm != NULL)
    {
      free(bm->map);
      free(bm);
    }
  potrace_param_free(param);
  return(st);
}

void rotate_point(int &x, int &y, int midX, int midY, double rotation)
{
// create 2D rotation matrix
  float sinval = sin(rotation); // no use of sincos()-function for compatibility, no performace bottleneck anymore anyway
  float cosval = cos(rotation);
  float m11 = cosval;
  float m12 = sinval;
  float m21 = -sinval;
  float m22 = cosval;

  int dX = x - midX;
  int dY = y - midY;

  int diffX = dX * m11 + dY * m21;
  int diffY = dX * m12 + dY * m22;

  x = midX+diffX;
  y = midY+diffY;
}

void rotate_coordinate_box(box_t &coordinate_box,double rotation,int width,int height)
{
  int midX = width/2;
  int midY = height/2;
  int x1 = coordinate_box.x1;
  int y1 = coordinate_box.y1;
  int x2 = coordinate_box.x2;
  int y2 = coordinate_box.y2;

  rotate_point(x1,y1,midX,midY,rotation);
  rotate_point(x2,y1,midX,midY,rotation);
  rotate_point(x1,y2,midX,midY,rotation);
  rotate_point(x2,y2,midX,midY,rotation);

  coordinate_box.x1 = std::min(x1,x2);
  coordinate_box.x2 = std::max(x1,x2);
  coordinate_box.y1 = std::min(y1,y2);
  coordinate_box.y2 = std::max(y1,y2);
}

void split_fragments_and_assemble_structure_record(
    std::vector<atom_t> &atom,
    int n_atom,
    std::vector<bond_t> &bond,
    int n_bond,
    const std::vector<box_t> &boxes,
    int l, int k,
    int resolution,
    int res_iter,
    const std::string &output_image_file_prefix,
    Image &image,
    Image &orig_box,
    int real_font_width, int real_font_height,
    double thickness,
    double avg_bond_length,
    const std::map<std::string, std::string> &superatom,
    int real_atoms, int real_bonds,
    int bond_max_type,
    double box_scale, double page_scale, double rotation, int unpaper_dx, int unpaper_dy,
    std::string output_format,
    const std::string &embedded_format,
    bool is_reaction,
    bool show_confidence,
    bool show_resolution_guess,
    bool show_page,
    bool show_coordinates,
    bool show_avg_bond_length,
    std::vector<std::vector<std::string> > &array_of_structures,
    std::vector<std::vector<double> > &array_of_avg_bonds,
    std::vector<std::vector<double> > &array_of_ind_conf,
    std::vector<std::vector<Image> > &array_of_images,
    std::vector<std::vector<box_t> > &array_of_boxes,
    int &total_boxes,
    double &total_confidence,
    int n_letters,
    bool show_learning,
    int resolution_iteration,
    bool verbose,
    const std::vector<bracket_t>&  brackets)
{
  std::vector<atom_t> frag_atom;
  std::vector<bond_t> frag_bond;

  if (real_atoms > MIN_A_COUNT && real_atoms < MAX_A_COUNT && real_bonds < MAX_B_COUNT && bond_max_type>0 && bond_max_type<5)
    {
      int num_frag;
      num_frag = resolve_bridge_bonds(atom, n_atom, bond, n_bond, 2 * thickness, avg_bond_length, superatom, verbose);
      collapse_bonds(atom, bond, n_bond, avg_bond_length / 4);
      collapse_atoms(atom, bond, n_atom, n_bond, 3);
      remove_zero_bonds(bond, n_bond, atom);
      extend_terminal_bond_to_bonds(atom, bond, n_bond, avg_bond_length, 7, 0);

      remove_small_terminal_bonds(bond, n_bond, atom, avg_bond_length);
      n_bond = reconnect_fragments(bond, n_bond, atom, avg_bond_length);

      collapse_atoms(atom, bond, n_atom, n_bond, 1);
      mark_terminal_atoms(bond, n_bond, atom, n_atom);
      const std::vector<std::vector<int> > &frags = find_fragments(bond, n_bond, atom);
      std::vector<fragment_t> fragments = populate_fragments(frags, atom);
      std::sort(fragments.begin(), fragments.end(), comp_fragments);
      for (unsigned int i = 0; i < fragments.size(); i++)
        {
          if (verbose)
            std::cout << "Considering fragment #" << i + 1 << " " << fragments[i].x1 << "x" << fragments[i].y1 << "-" << fragments[i].x2 << "x"
                      << fragments[i].y2 << ", atoms: " << fragments[i].atom.size() << '.' << std::endl;

          if (fragments[i].atom.size() > MIN_A_COUNT)
            {
              frag_atom.clear();
              for (int a = 0; a < n_atom; a++)
                {
                  frag_atom.push_back(atom[a]);
                  frag_atom[a].exists = false;
                }

              for (unsigned int j = 0; j < fragments[i].atom.size(); j++)
                frag_atom[fragments[i].atom[j]].exists = atom[fragments[i].atom[j]].exists;

              frag_bond.clear();
              for (int b = 0; b < n_bond; b++)
                {
                  frag_bond.push_back(bond[b]);
                }

              remove_zero_bonds(frag_bond, n_bond, frag_atom);

              double confidence = 0;
              molecule_statistics_t molecule_statistics;
              int page_number = l + 1;
              box_t coordinate_box,rel_box;
	       if (fragments.size()>1)
		{
		  coordinate_box.x1 = (int) (-(double)page_scale * unpaper_dx + (double) page_scale * boxes[k].x1 + (double) page_scale * box_scale * fragments[i].x1 - (double) page_scale * FRAME);
		  coordinate_box.y1 = (int) (-(double)page_scale * unpaper_dy + (double) page_scale * boxes[k].y1 + (double) page_scale * box_scale * fragments[i].y1 - (double) page_scale * FRAME);
		  coordinate_box.x2 = (int) (-(double)page_scale * unpaper_dx + (double) page_scale * boxes[k].x1 + (double) page_scale * box_scale * fragments[i].x2 - (double) page_scale * FRAME);
		  coordinate_box.y2 = (int) (-(double)page_scale * unpaper_dy + (double) page_scale * boxes[k].y1 + (double) page_scale * box_scale * fragments[i].y2 - (double) page_scale * FRAME);
		  //rotate_coordinate_box(coordinate_box,rotation,image.columns(),image.rows());
		  rel_box.x1 = (int)((double)boxes[k].x1 + (double) box_scale * fragments[i].x1 - FRAME);
		  rel_box.y1 = (int)((double)boxes[k].y1 + (double) box_scale * fragments[i].y1 - FRAME);
		  rel_box.x2 = (int)((double)boxes[k].x1 + (double) box_scale * fragments[i].x2 - FRAME);
		  rel_box.y2 = (int)((double)boxes[k].y1 + (double) box_scale * fragments[i].y2 - FRAME);
		}
		else
		{
		  coordinate_box.x1 = (int) (-(double)page_scale * unpaper_dx + (double) page_scale * boxes[k].x1);
		  coordinate_box.y1 = (int) (-(double)page_scale * unpaper_dy + (double) page_scale * boxes[k].y1);
		  coordinate_box.x2 = (int) (-(double)page_scale * unpaper_dx + (double) page_scale * boxes[k].x2);
		  coordinate_box.y2 = (int) (-(double)page_scale * unpaper_dy + (double) page_scale * boxes[k].y2);
		  //rotate_coordinate_box(coordinate_box,rotation,image.columns(),image.rows());
		  rel_box.x1 = boxes[k].x1;
		  rel_box.y1 = boxes[k].y1;
		  rel_box.x2 = boxes[k].x2;
		  rel_box.y2 = boxes[k].y2;
		  }

              if (verbose)
                std::cout << "Coordinate box: " << coordinate_box.x1 << "x" << coordinate_box.y1 << "-" << coordinate_box.x2 << "x"
                          << coordinate_box.y2 << "." << std::endl;
	      if (is_reaction)
		output_format = SUBSTITUTE_REACTION_FORMAT;

              std::string structure =
                get_formatted_structure(frag_atom, frag_bond, n_bond, output_format, embedded_format,
                                        molecule_statistics, confidence,
                                        show_confidence, avg_bond_length, page_scale * box_scale * avg_bond_length,
                                        show_avg_bond_length,
                                        show_resolution_guess ? &resolution : NULL,
                                        show_page ? &page_number : NULL,
                                        show_coordinates ? &coordinate_box : NULL, superatom, n_letters, show_learning, resolution_iteration, verbose,
					brackets);

              if (molecule_statistics.fragments > 0 && molecule_statistics.fragments < MAX_FRAGMENTS
		  && molecule_statistics.num_atoms>MIN_A_COUNT && molecule_statistics.num_bonds>0
		  )
                {
		  if ((molecule_statistics.rings56 > 0 || molecule_statistics.num_organic_non_carbon_atoms > 0)
		      && molecule_statistics.num_bonds>MIN_B_COUNT
		      && molecule_statistics.num_small_angles < 3
		      && avg_bond_length > real_font_height )
		    {
		      array_of_structures[res_iter].push_back(structure);
		      array_of_avg_bonds[res_iter].push_back(page_scale * box_scale * avg_bond_length);
		      array_of_ind_conf[res_iter].push_back(confidence);
		      array_of_boxes[res_iter].push_back(rel_box);

		      if (!output_image_file_prefix.empty())
			{
			  Image tmp = image;
			  if (!is_reaction)
			    {
			      Geometry geometry =
				(fragments.size() > 1) ? Geometry(box_scale * fragments[i].x2 - box_scale * fragments[i].x1, //
								  box_scale * fragments[i].y2 - box_scale * fragments[i].y1, //
								  boxes[k].x1 + box_scale * fragments[i].x1 - FRAME , //
								  boxes[k].y1 + box_scale * fragments[i].y1 - FRAME )
				: Geometry(boxes[k].x2 - boxes[k].x1, boxes[k].y2 - boxes[k].y1, boxes[k].x1, boxes[k].y1);

			      try
				{
				  tmp.crop(geometry);
				}
			      catch (...)
				{
				  tmp = orig_box;
				}
			    }
			  array_of_images[res_iter].push_back(tmp);
			}
		    }
                  total_boxes++;
                  total_confidence += confidence;
		  if (verbose)
                    std::cout << "Result: " << res_iter << " " << structure << " " << confidence << std::endl;
                }
            }
        }
    }
}

int count_recognized_chars(std::vector<atom_t>  &atom, std::vector<bond_t>& bond)
{
  std::string char_filter = "oOcCNHsSBMeEXYZRPp23456789AF";
  std::set<int> atoms;
  for (int i=0; i<bond.size(); i++)
    if (bond[i].exists)
      {
	atoms.insert(bond[i].a);
	atoms.insert(bond[i].b);
      }
  int r = 0;
  for (std::set<int>::iterator a = atoms.begin(); a != atoms.end(); a++)
    for (int i = 0; i < atom[*a].label.size(); i++)
      {
        if (char_filter.find(atom[*a].label[i]) != std::string::npos)
	  r++;
      }
  return r;
}

Image process_pdf_page(poppler::document* doc,  poppler::page_renderer &r, int l, int resolution)
{
  poppler::page* p = doc->create_page(l);
  poppler::image im = r.render_page(p, resolution, resolution);
  Image image(Geometry(im.width(), im.height()), "white");
  image.modifyImage();
  image.type(TrueColorType);
  const char *d = im.const_data();
  int bytes_per_row = im.bytes_per_row();
  int bytes_per_pixel = bytes_per_row / im.width();
  for (int row = 0; row < im.height(); ++row)
    {
      for (int col = 0; col < im.width(); ++col)
	{
	  int offset = row * bytes_per_row + col * bytes_per_pixel;
	  int r = -1, g = -1, b = -1, a = -1;
	      switch(im.format())
		{
		case poppler::image::format_mono : r = g = b = int(d[offset]); break;
		case poppler::image::format_rgb24  : r = int(d[offset]); g = int(d[offset+1]); b = int(d[offset+2]); break;
		case poppler::image::format_argb32 : r = int(d[offset]); g = int(d[offset+1]); b = int(d[offset+2]); a = int(d[offset+3]); break;
		}
	      if (r >= 0 && g >= 0 && b >= 0)
		image.pixelColor(col, row, ColorRGB(double(r) / 128, double(g) / 128, double(b) / 128));
	}
    }
  return image;
}

extern job_t *OCR_JOB;
extern job_t *JOB;

#ifdef OSRA_LIB
int global_init_state;
#endif

// Function: osra_init()
//
// Initialises OSRA library. Should be called at e.g. program startup. This function is automatically called for both SO library and CLI utility.
// See this section for details about library init/cleanup: http://www.faqs.org/docs/Linux-HOWTO/Program-Library-HOWTO.html#INIT-AND-CLEANUP
// Below attribute marker is GNU compiler specific.
void __attribute__ ((constructor)) osra_init()
{
  // Necessary for GraphicsMagick-1.3.8 according to http://www.graphicsmagick.org/1.3/NEWS.html#january-21-2010:
  MagickLib::InitializeMagick(NULL);

  osra_ocr_init();

#ifdef OSRA_LIB
  global_init_state = osra_openbabel_init();

  if (global_init_state != 0)
    std::cerr << "OpenBabel initialization failure." << std::endl;
#endif

  srand(1);
}

// Function: osra_destroy()
//
// Releases all resources allocated by OSRA library. Should be called at e.g. program exit. This function is automatically called for both SO library and CLI utility.
// See this section for details about library init/cleanup: http://www.faqs.org/docs/Linux-HOWTO/Program-Library-HOWTO.html#INIT-AND-CLEANUP
// Below attribute marker is GNU compiler specific.
void __attribute__ ((destructor)) osra_destroy()
{
#ifdef OSRA_LIB
  MagickLib::DestroyMagick();
#endif
  osra_ocr_destroy();
}

int osra_process_image(
#ifdef OSRA_LIB
  const char *image_data,
  int image_length,
  std::ostream &structure_output_stream,
#else
  const std::string &input_file,
  const std::string &output_file,
#endif
  int rotate,
  bool invert,
  int input_resolution,
  double threshold,
  int do_unpaper,
  bool jaggy,
  bool adaptive_option,
  std::string output_format,
  std::string embedded_format,
  bool show_confidence,
  bool show_resolution_guess,
  bool show_page,
  bool show_coordinates,
  bool show_avg_bond_length,
  bool show_learning,
  const std::string &osra_dir,
  const std::string &spelling_file,
  const std::string &superatom_file,
  bool debug,
  bool verbose,
  const std::string &output_image_file_prefix,
  const std::string &resize,
  const std::string &preview
)
{
#ifdef OSRA_LIB
  if (global_init_state != 0) return global_init_state;
#endif

  std::transform(output_format.begin(), output_format.end(), output_format.begin(), ::tolower);
  std::transform(embedded_format.begin(), embedded_format.end(), embedded_format.begin(), ::tolower);

  std::map<std::string, std::string> spelling, superatom;
  int err = load_superatom_spelling_maps(spelling, superatom, osra_dir, spelling_file, superatom_file, verbose);
  if (err != 0) return err;

  std::string type;

#ifdef OSRA_LIB
  Blob blob(image_data, image_length);
#endif

  try
    {
      Image image_typer;
#ifdef OSRA_LIB
      image_typer.ping(blob);
#else
      image_typer.ping(input_file);
#endif
      type = image_typer.magick();
    }
  catch (...)
    {
      // Unfortunately, GraphicsMagick does not throw exceptions in all cases, so it behaves inconsistent, see
      // https://sourceforge.net/tracker/?func=detail&aid=3022955&group_id=40728&atid=428740
    }

  //int stderr_copy = dup(2);
  //fclose(stderr);

  int page = 1;
  poppler::document* poppler_doc = NULL;
  if (type.empty() || type == "PDF" || type == "PS")
    {
#ifdef OSRA_LIB
      poppler_doc = poppler::document::load_from_raw_data(image_data, image_length);
#else
      poppler_doc = poppler::document::load_from_file(input_file);
#endif
    }
  if (poppler_doc)
    {
      page = poppler_doc->pages();
      type = "PDF";
    }
  else if (type == "PDF" || type == "PS")
    {
      type.clear();
    }
  else if (!type.empty() && type != "PDF")
    {
#ifdef OSRA_LIB
      page = count_pages(blob);
#else
      page = count_pages(input_file);
#endif
    }
  // dup2(stderr_copy, 2);
  //close(stderr_copy);

  if (type.empty())
    {
#ifdef OSRA_LIB
      std::cerr << "Cannot detect blob image type" << std::endl;
#else
      std::cerr << "Cannot open file \"" << input_file << '"' << std::endl;
#endif
      return ERROR_UNKNOWN_IMAGE_TYPE;
    }

  if (verbose)
    std::cout << "Image type: " << type << '.' << std::endl;

#ifndef OSRA_LIB
  std::ofstream outfile;

  if (!output_file.empty())
    {
      outfile.open(output_file.c_str(), std::ios::out | std::ios::trunc);
      if (outfile.bad() || !outfile.is_open())
        {
          std::cerr << "Cannot open file \"" << output_file << "\" for output" << std::endl;
          return ERROR_OUTPUT_FILE_OPEN_FAILED;
        }
    }
#endif


  if (show_coordinates && rotate != 0)
    {
      std::cerr << "Showing the box coordinates is currently not supported together with image rotation and is therefore disabled." << std::endl;
#ifdef OSRA_LIB
      return ERROR_ILLEGAL_ARGUMENT_COMBINATION;
#else
      show_coordinates = false;
#endif
    }

  if (!embedded_format.empty() && !(output_format == "sdf" && (embedded_format == "inchi" || embedded_format == "smi"
                                    || embedded_format == "can")))
    {
      std::cerr << "Embedded format option is only possible if output format is SDF and option can have only inchi, smi, or can values." << std::endl;
      return ERROR_ILLEGAL_ARGUMENT_COMBINATION;
    }

  // This will hide the output "Warning: non-positive median line gap" from GOCR. Remove after this is fixed:
  fclose(stderr);
  OpenBabel::obErrorLog.StopLogging();

  bool is_reaction = false;
  if (output_format == "cmlr" || output_format == "rsmi" || output_format =="rxn")
    is_reaction = true;


  std::vector<std::vector<std::string> > pages_of_structures(page, std::vector<std::string> (0));
  std::vector<std::vector<Image> > pages_of_images(page, std::vector<Image> (0));
  std::vector<std::vector<double> > pages_of_avg_bonds(page, std::vector<double> (0));
  std::vector<std::vector<double> > pages_of_ind_conf(page, std::vector<double> (0));
  std::vector<std::vector<box_t> > pages_of_boxes(page, std::vector<box_t> (0));
  std::vector<std::vector<arrow_t> > arrows(page, std::vector<arrow_t>(0));
  std::vector<std::vector<plus_t> > pluses(page, std::vector<plus_t>(0));

  int total_structure_count = 0;
  int num_resolutions = NUM_RESOLUTIONS;
  if (input_resolution != 0)
    num_resolutions = 1;
  std::vector<double> array_of_confidence(num_resolutions, 0);
  std::vector<int> boxes_per_res(num_resolutions,0);
  std::vector<int> select_resolution(num_resolutions, input_resolution);
  set_select_resolution(select_resolution,input_resolution);
  std::vector<std::vector<std::vector<std::string> > > array_of_structures_page(
      page, std::vector<std::vector<std::string> >(num_resolutions));
  std::vector<std::vector<std::vector<double> > > array_of_avg_bonds_page(
      page, std::vector<std::vector<double> >(num_resolutions));
  std::vector<std::vector<std::vector<double> > > array_of_ind_conf_page(
      page, std::vector<std::vector<double> >(num_resolutions));
  std::vector<std::vector<std::vector<Image> > > array_of_images_page(
      page, std::vector<std::vector<Image> > (num_resolutions));
  std::vector<std::vector<std::vector<box_t> > > array_of_boxes_page(
      page, std::vector<std::vector<box_t> >(num_resolutions));

//#pragma omp parallel for default(shared) private(OCR_JOB,JOB)
  for (int l = 0; l < page; l++)
    {
      Image image;
      double page_scale=1;
      poppler::page_renderer poppler_renderer;

      int ttt = 0;

      if (verbose)
        std::cout << "Processing page " << (l+1) << " out of " << page << "..." << std::endl;

      if (type == "PDF" || type == "PS")
        page_scale *= (double) 72 / input_resolution;


      if (poppler_doc) // process PDF and PS files
	{
	  int resolution = input_resolution;
	  if (resolution == 0)
	    resolution = 300;
	  image = process_pdf_page(poppler_doc, poppler_renderer, l, resolution);
	}
      else
	{
#ifdef OSRA_LIB
	  image.read(blob);
#else
          std::ostringstream pname;
	  pname << input_file << "[" << l << "]";
#pragma omp critical
	  {
	    image.read(pname.str());
	  }
#endif
	}
      if (l == 0 && !preview.empty())
	{
	  try
	    {
	      image.write(preview);
	    }
	  catch(...)
	    {}
	}

      image.modifyImage();
      bool adaptive = convert_to_gray(image, invert, adaptive_option, verbose);

      std::vector<std::vector<std::string> > array_of_structures(num_resolutions);
      std::vector<std::vector<double> > array_of_avg_bonds(num_resolutions), array_of_ind_conf(num_resolutions);
      std::vector<std::vector<Image> > array_of_images(num_resolutions);
      std::vector<std::vector<box_t> > array_of_boxes(num_resolutions);

      if (input_resolution > 300)
        {
          int percent = (100 * 300) / input_resolution;
          std::ostringstream scale;
          scale << percent << "%";
          image.scale(scale.str());
          page_scale /= (double) percent / 100;
        }

      if (verbose)
        std::cout << "Input resolutions are " << select_resolution << std::endl;

      ColorGray bgColor = getBgColor(image);
      if (rotate != 0)
        {
          image.backgroundColor(bgColor);
          image.rotate(rotate);
        }

      double rotation = 0;
      int unpaper_dx = 0;
      int unpaper_dy = 0;
      for (int i = 0; i < do_unpaper; i++)
        {
          double radians=0;
          int dx=0, dy=0;
          unpaper(image,radians,dx,dy);
          rotation +=radians;
          unpaper_dx +=dx;
          unpaper_dy +=dy;
        }

      // 0.1 is used for THRESHOLD_BOND here to allow for farther processing.
      std::list<std::list<std::list<point_t> > > clusters = find_segments(image, 0.1, bgColor, adaptive, is_reaction, arrows[l], pluses[l], verbose);

      if (verbose)
        std::cout << "Number of clusters: " << clusters.size() << '.' << std::endl;

      std::vector<box_t> boxes;
      std::set<std::pair<int, int> > brackets;
      int n_boxes = prune_clusters(clusters, boxes, brackets);
      std::sort(boxes.begin(), boxes.end(), comp_boxes);

      if (verbose)
        std::cout << "Number of boxes: " << boxes.size() << '.' << std::endl;


      for (int res_iter = 0; res_iter < num_resolutions; res_iter++)
        {
          int total_boxes = 0;
          double total_confidence = 0;

          int resolution = select_resolution[res_iter];
          int working_resolution = resolution;
          if (resolution > 300)
            working_resolution = 300;

          double THRESHOLD_BOND = set_threshold(threshold,resolution);

          int max_font_height = MAX_FONT_HEIGHT * working_resolution / 150;
          int max_font_width = MAX_FONT_WIDTH * working_resolution / 150;
          bool thick = true;
          if (resolution < 150)
            thick = false;
          else if (resolution == 150 && !jaggy)
            thick = false;

          //Image dbg = image;
          //dbg.modifyImage();
          //dbg.backgroundColor("white");
          //dbg.erase();
          //dbg.type(TrueColorType);
          for (int k = 0; k < n_boxes; k++)
            if ((boxes[k].x2 - boxes[k].x1) > max_font_width && (boxes[k].y2 - boxes[k].y1) > max_font_height
                && !boxes[k].c.empty() && ((boxes[k].x2 - boxes[k].x1) > 2 * max_font_width || (boxes[k].y2
                                           - boxes[k].y1) > 2 * max_font_height))
              {
                int n_atom = 0, n_bond = 0, n_letters = 0, n_label = 0;
                std::vector<atom_t> atom;
                std::vector<bond_t> bond;
                std::vector<letters_t> letters;
                std::vector<label_t> label;
                double box_scale = 1;
                Image orig_box(Geometry(boxes[k].x2 - boxes[k].x1 + 2 * FRAME, boxes[k].y2 - boxes[k].y1 + 2
                                        * FRAME), bgColor);

                for (unsigned int p = 0; p < boxes[k].c.size(); p++)
                  {
                    int x = boxes[k].c[p].x;
                    int y = boxes[k].c[p].y;
                    ColorGray color = image.pixelColor(x, y);
                    //dbg.pixelColor(x, y, color);
                    orig_box.pixelColor(x - boxes[k].x1 + FRAME, y - boxes[k].y1 + FRAME, color);
                  }


                int width = orig_box.columns();
                int height = orig_box.rows();
                Image thick_box;
                create_thick_box(orig_box,thick_box,width,height,resolution,working_resolution,box_scale,bgColor,THRESHOLD_BOND,res_iter,thick,jaggy);

                if (verbose)
                  std::cout << "Analysing box " << boxes[k].x1 << "x" << boxes[k].y1 << "-" << boxes[k].x2 << "x" << boxes[k].y2 << " using working resolution " << working_resolution << '.' << std::endl;

                Image box;
                if (thick)
                  box = thin_image(thick_box, THRESHOLD_BOND, bgColor);
                else
                  box = thick_box;
                potrace_state_t * const  st = raster_to_vector(box,bgColor,THRESHOLD_BOND,width,height,working_resolution);
                potrace_path_t const * const p = st->plist;
                n_atom = find_atoms(p, atom, bond, &n_bond,width,height);

                int real_font_width, real_font_height;
                n_letters = find_chars(p, orig_box, letters, atom, bond, n_atom, n_bond, height, width, bgColor,
                                       THRESHOLD_BOND, max_font_width, max_font_height, real_font_width, real_font_height,verbose);
                if (verbose)
                  std::cout << "Number of atoms: " << n_atom << ", bonds: " << n_bond << ", " << n_letters << " letters: " << n_letters << " " << letters << " after find_atoms()" << std::endl;

                double avg_bond_length = percentile75(bond, n_bond, atom);

                double max_area = avg_bond_length * 5;
                if (thick)
                  max_area = avg_bond_length;
                n_letters = find_plus_minus(p, orig_box, bgColor, THRESHOLD_BOND, letters, atom, bond, n_atom, n_bond, height, width,
                                            real_font_height, real_font_width, n_letters, avg_bond_length);
                n_atom = find_small_bonds(p, atom, bond, n_atom, &n_bond, max_area, avg_bond_length / 2, 5);

		//remove_small_bonds_in_chars(atom,bond,letters);

                find_old_aromatic_bonds(p, bond, n_bond, atom, n_atom, avg_bond_length);

                if (verbose)
                  std::cout << "Number of atoms: " << n_atom << ", bonds: " << n_bond << ", " << n_letters << "letters: " << letters << " after find_old_aromatic_bonds()" << std::endl;

                double dist = 3.;
                if (working_resolution < 150)
                  dist = 2;

                double thickness = skeletize(atom, bond, n_bond, box, THRESHOLD_BOND, bgColor, dist, avg_bond_length);
                remove_disconnected_atoms(atom, bond, n_atom, n_bond);
                collapse_atoms(atom, bond, n_atom, n_bond, 3);
                remove_zero_bonds(bond, n_bond, atom);

		n_bond = find_wavy_bonds(bond,n_bond,atom,avg_bond_length);
		//				if (ttt++ == 0)  debug_image(orig_box, atom, n_atom, bond, n_bond, "tmp.png");
                n_letters = find_fused_chars(bond, n_bond, atom, letters, n_letters, real_font_height,
                                             real_font_width, 0, orig_box, bgColor, THRESHOLD_BOND, 3, verbose);

                n_letters = find_fused_chars(bond, n_bond, atom, letters, n_letters, real_font_height,
                                             real_font_width, '*', orig_box, bgColor, THRESHOLD_BOND, 5, verbose);

                flatten_bonds(bond, n_bond, atom, 3);
                remove_zero_bonds(bond, n_bond, atom);
                avg_bond_length = percentile75(bond, n_bond, atom);

                if (verbose)
                  std::cout << "Average bond length: " << avg_bond_length << std::endl;

                double max_dist_double_bond = dist_double_bonds(atom, bond, n_bond, avg_bond_length);
                n_bond = double_triple_bonds(atom, bond, n_bond, avg_bond_length, n_atom, max_dist_double_bond);
                n_atom = find_dashed_bonds(p, atom, bond, n_atom, &n_bond, std::max(MAX_DASH, int(avg_bond_length / 3)),
                                           avg_bond_length, orig_box, bgColor, THRESHOLD_BOND, thick, avg_bond_length, letters);

                n_letters = remove_small_bonds(bond, n_bond, atom, letters, n_letters, real_font_height,
                                               MIN_FONT_HEIGHT, avg_bond_length);

		n_letters = find_numbers(p, orig_box, letters, atom, bond, n_atom, n_bond, height, width, bgColor,
					 THRESHOLD_BOND, n_letters);

                dist = 4.;
                if (working_resolution < 300)
                  dist = 3;
                if (working_resolution < 150)
                  dist = 2;

                n_bond = fix_one_sided_bonds(bond, n_bond, atom, dist, avg_bond_length);

                n_letters = clean_unrecognized_characters(bond, n_bond, atom, real_font_height, real_font_width, 4,
		          letters, n_letters);

                thickness = find_wedge_bonds(thick_box, atom, n_atom, bond, n_bond, bgColor, THRESHOLD_BOND,
                                             max_dist_double_bond, avg_bond_length, 3, 1);

                n_label = assemble_labels(letters, n_letters, label);

                if (verbose)
                  std::cout << n_label << " labels: " << label << " after assemble_labels()" << std::endl;

                remove_disconnected_atoms(atom, bond, n_atom, n_bond);

                collapse_atoms(atom, bond, n_atom, n_bond, thickness);

                remove_zero_bonds(bond, n_bond, atom);

                flatten_bonds(bond, n_bond, atom, 2 * thickness);

                remove_zero_bonds(bond, n_bond, atom);

                avg_bond_length = percentile75(bond, n_bond, atom);

                collapse_double_bonds(bond, n_bond, atom, max_dist_double_bond);

                extend_terminal_bond_to_label(atom, letters, n_letters, bond, n_bond, label, n_label, avg_bond_length / 2,
					      thickness, max_dist_double_bond);

                remove_disconnected_atoms(atom, bond, n_atom, n_bond);
                collapse_atoms(atom, bond, n_atom, n_bond, thickness);
                collapse_doubleup_bonds(bond, n_bond);

                remove_zero_bonds(bond, n_bond, atom);
                flatten_bonds(bond, n_bond, atom, thickness);
                remove_zero_bonds(bond, n_bond, atom);
                remove_disconnected_atoms(atom, bond, n_atom, n_bond);

                extend_terminal_bond_to_bonds(atom, bond, n_bond, avg_bond_length, 2 * thickness, max_dist_double_bond);

                std::vector<bracket_t> bracket_boxes;
		remove_bracket_atoms(atom, n_atom, bond, n_bond, brackets, thickness, boxes[k].x1, boxes[k].y1, box_scale, real_font_width, real_font_height, bracket_boxes);
		remove_zero_bonds(bond, n_bond, atom);
		remove_vertical_bonds_close_to_brackets(bracket_boxes, atom, bond, n_bond, thickness, avg_bond_length);
		remove_zero_bonds(bond, n_bond, atom);
		flatten_bonds(bond, n_bond, atom, 2*thickness);
		assign_labels_to_brackets(bracket_boxes, label, n_label, letters, n_letters, real_font_width, real_font_height);

                collapse_atoms(atom, bond, n_atom, n_bond, 3);

                remove_zero_bonds(bond, n_bond, atom);
                flatten_bonds(bond, n_bond, atom, 5);
                remove_zero_bonds(bond, n_bond, atom);

                n_letters = clean_unrecognized_characters(bond, n_bond, atom, real_font_height, real_font_width, 0,
                            letters, n_letters);
		int recognized_chars = count_recognized_chars(atom,bond);


                assign_charge(atom, bond, n_atom, n_bond, spelling, superatom, debug);
                find_up_down_bonds(bond, n_bond, atom, thickness);
                int real_atoms = count_atoms(atom, n_atom);
                int bond_max_type = 0;
                int real_bonds = count_bonds(bond, n_bond,bond_max_type);

                if (verbose)
                  std::cout << "Final number of atoms: " << real_atoms << ", bonds: " << real_bonds << ", chars: " << n_letters << '.' << std::endl;


                split_fragments_and_assemble_structure_record(atom,n_atom,bond,n_bond,boxes,
							      l,k,resolution,res_iter,output_image_file_prefix,image,orig_box,real_font_width,real_font_height,
							      thickness,avg_bond_length,superatom,real_atoms,real_bonds,bond_max_type,
							      box_scale,page_scale,rotation,unpaper_dx,unpaper_dy,output_format,embedded_format,is_reaction,show_confidence,
							      show_resolution_guess,show_page,show_coordinates, show_avg_bond_length,array_of_structures,
							      array_of_avg_bonds,array_of_ind_conf,array_of_images,array_of_boxes,total_boxes,total_confidence,
							      recognized_chars,show_learning,res_iter,verbose, bracket_boxes);

                if (st != NULL)
                  potrace_state_free(st);
              }
	  array_of_confidence[res_iter] += total_confidence;
	  boxes_per_res[res_iter] += total_boxes;
          //dbg.write("debug.png");
        }

      #pragma omp critical
      {
         if (show_learning)
	  for (int j = 0; j < num_resolutions; j++)
	    for (unsigned int i = 0; i < array_of_structures[j].size(); i++)
	    {
	      pages_of_structures[l].push_back(array_of_structures[j][i]);
	      if (!output_image_file_prefix.empty())
		pages_of_images[l].push_back(array_of_images[j][i]);
	      pages_of_avg_bonds[l].push_back(array_of_avg_bonds[j][i]);
	      pages_of_ind_conf[l].push_back(array_of_ind_conf[j][i]);
	      pages_of_boxes[l].push_back(array_of_boxes[j][i]);
	      total_structure_count++;
	    }
	 else
	   for (int j = 0; j < num_resolutions; j++)
	     {
                array_of_structures_page[l][j] = array_of_structures[j];
		if (!output_image_file_prefix.empty())
		  array_of_images_page[l][j] = array_of_images[j];
		array_of_avg_bonds_page[l][j] = array_of_avg_bonds[j];
		array_of_ind_conf_page[l][j] = array_of_ind_conf[j];
		array_of_boxes_page[l][j] = array_of_boxes[j];
	      }

       }
     }

    double max_conf = -FLT_MAX;
    int max_res = 0;
    for (int i = 0; i < num_resolutions; i++)
        {
          if (boxes_per_res[i] > 0 && array_of_confidence[i]/boxes_per_res[i] > max_conf)
            {
              max_conf = array_of_confidence[i]/boxes_per_res[i];
              max_res = i;
            }
        }
      for (int i = 0; i < num_resolutions; i++)
	if (boxes_per_res[i] > 0 && array_of_confidence[i]/boxes_per_res[i] == max_conf && select_resolution[i] == 300) // second 300 dpi is without thinning
	  {
	    max_res = i;
	    break;
	  }

	if (!show_learning)
         for (int l = 0; l < page; l++)
	    {
	      pages_of_structures[l] = array_of_structures_page[l][max_res];
	      if (!output_image_file_prefix.empty())
		pages_of_images[l] = array_of_images_page[l][max_res];
	      pages_of_avg_bonds[l] = array_of_avg_bonds_page[l][max_res];
	      pages_of_ind_conf[l] = array_of_ind_conf_page[l][max_res];
	      pages_of_boxes[l] = array_of_boxes_page[l][max_res];
	      total_structure_count += array_of_structures_page[l][max_res].size();
	    }

  double best_bond = 0;

  //if (total_structure_count >= STRUCTURE_COUNT)
  //  find_limits_on_avg_bond(best_bond, pages_of_avg_bonds, pages_of_ind_conf);

  // If multiple pages are processed at several  resolutions different pages
  // may be processed at different resolutions leading to a seemingly different average bond length
  // Currently multi-page documents (PDF and PS) are all processed at the same resolution
  // and single-page images have all structures on the page at the same resolution

  //cout << min_bond << " " << max_bond << endl;

#ifdef OSRA_LIB
  std::ostream &out_stream = structure_output_stream;
#else
  std::ostream &out_stream = outfile.is_open() ? outfile : std::cout;
#endif

  // For Andriod version we will find the structure with maximum confidence value, as the common usecase for Andriod is to analyse the
  // image (taken by embedded photo camera) that usually contains just one molecule:
  double max_confidence = -FLT_MAX;
  int l_index = 0;
  int i_index = 0;
  int image_count = 0;

  for (int l = 0; l < page; l++)
    {
      for (unsigned int i = 0; i < pages_of_structures[l].size(); i++)
	if (best_bond == 0 || (pages_of_avg_bonds[l][i] > best_bond/2 && pages_of_avg_bonds[l][i] < 2*best_bond))
	  {
	    if (pages_of_ind_conf[l][i] > max_confidence)
	      {
		max_confidence = pages_of_ind_conf[l][i];
		l_index = l;
		i_index = i;
	      }

	    if (output_format != "mol" && !is_reaction)
	      {
		out_stream << pages_of_structures[l][i];

		// Dump this structure into a separate file:
		if (!output_image_file_prefix.empty())
		  {
                    std::ostringstream fname;
		    fname << output_image_file_prefix << image_count << ".png";
		    image_count++;
		    if (fname.str() != "")
		      {
			Image tmp = pages_of_images[l][i];
			if (resize != "")
			  {
			    tmp.scale(resize);
			  }
			tmp.write(fname.str());
		      }
		  }
	      }
	  }
      if (is_reaction && !arrows[l].empty())
	 {
           std::vector<std::string> reactions;
           std::vector<box_t> rbox;
	   arrange_reactions(arrows[l], pages_of_boxes[l], pluses[l], reactions, rbox, pages_of_structures[l],output_format);
	   for (int k=0; k<reactions.size(); k++)
	     {
               out_stream << reactions[k] << std::endl;

	       if (!output_image_file_prefix.empty())
		 {
                   std::ostringstream fname;
		   fname << output_image_file_prefix << image_count << ".png";
		   image_count++;
		   if (fname.str() != "")
		     {
		       Image tmp = pages_of_images[l][k];
		       Geometry geometry = Geometry(rbox[k].x2 - rbox[k].x1, rbox[k].y2 - rbox[k].y1, rbox[k].x1, rbox[k].y1);
		       tmp.crop(geometry);
		       if (resize != "")
			 {
			   tmp.scale(resize);
			 }
		       tmp.write(fname.str());
		     }
		 }
	     }
	 }
    }
  // Output the structure with maximum confidence value:
  if (output_format == "mol")
    {
      out_stream << pages_of_structures[l_index][i_index];
      if (!output_image_file_prefix.empty())
	{
          std::ostringstream fname;
	  fname << output_image_file_prefix  << ".png";
	  if (fname.str() != "")
	    {
	      Image tmp = pages_of_images[l_index][i_index];
	      if (resize != "")
		{
		  tmp.scale(resize);
		}
	      tmp.write(fname.str());
	    }
	}
    }

  out_stream.flush();

#ifndef OSRA_LIB
  if (!output_file.empty())
    outfile.close();
#endif

  return 0;
}

int hack_osra_process_image(
  const char *image_data,
  int image_length,
  std::vector<std::string> rgroup_vars,
  std::vector<std::string> rgroup_values,
  const std::string &input_file,
  const std::string &output_file,
  int rotate,
  bool invert,
  int input_resolution,
  double threshold,
  int do_unpaper,
  bool jaggy,
  bool adaptive_option,
  std::string output_format,
  std::string embedded_format,
  bool show_confidence,
  bool show_resolution_guess,
  bool show_page,
  bool show_coordinates,
  bool show_avg_bond_length,
  bool show_learning,
  const std::string &osra_dir,
  const std::string &spelling_file,
  const std::string &superatom_file,
  bool debug,
  bool verbose,
  const std::string &output_image_file_prefix,
  const std::string &resize,
  const std::string &preview
//  std::vector<std::vector<std::string, std::string>> rgroup
)
{
bool from_file = false;
//char image_data = 'a';
//int image_length = 4;
int global_init_state;


if (from_file == true){
  if (global_init_state != 0) return global_init_state;
}

  std::transform(output_format.begin(), output_format.end(), output_format.begin(), ::tolower);
  std::transform(embedded_format.begin(), embedded_format.end(), embedded_format.begin(), ::tolower);

  std::map<std::string, std::string> spelling, superatom;
  int err = load_superatom_spelling_maps(spelling, superatom, osra_dir, spelling_file, superatom_file, verbose);
  if (err != 0) return err;

  std::string type;

if (from_file == true){
  Blob blob(image_data, image_length);
}

  try
    {
      Image image_typer;

      image_typer.ping(input_file);

      type = image_typer.magick();
    }
  catch (...)
    {
      // Unfortunately, GraphicsMagick does not throw exceptions in all cases, so it behaves inconsistent, see
      // https://sourceforge.net/tracker/?func=detail&aid=3022955&group_id=40728&atid=428740
    }

  //int stderr_copy = dup(2);
  //fclose(stderr);

  int page = 1;
  poppler::document* poppler_doc = NULL;
  if (type.empty() || type == "PDF" || type == "PS")
    {
if (from_file == true){
      poppler_doc = poppler::document::load_from_raw_data(image_data, image_length);
} else{
      poppler_doc = poppler::document::load_from_file(input_file);
}
    }
  if (poppler_doc)
    {
      page = poppler_doc->pages();
      type = "PDF";
    }
  else if (type == "PDF" || type == "PS")
    {
      type.clear();
    }
  else if (!type.empty() && type != "PDF")
    {

      page = count_pages(input_file);

    }
  // dup2(stderr_copy, 2);
  //close(stderr_copy);

  if (type.empty())
    {
if (from_file = true){
      std::cerr << "Cannot detect blob image type" << std::endl;
} else {
      std::cerr << "Cannot open file \"" << input_file << '"' << std::endl;
}
      return ERROR_UNKNOWN_IMAGE_TYPE;
    }

  if (verbose)
    std::cout << "Image type: " << type << '.' << std::endl;

std::ofstream outfile;

if (!output_file.empty())
{
  outfile.open(output_file.c_str(), std::ios::out | std::ios::trunc);
  if (outfile.bad() || !outfile.is_open())
    {
      std::cerr << "Cannot open file \"" << output_file << "\" for output" << std::endl;
      return ERROR_OUTPUT_FILE_OPEN_FAILED;
    }
}



  if (show_coordinates && rotate != 0)
    {
      std::cerr << "Showing the box coordinates is currently not supported together with image rotation and is therefore disabled." << std::endl;
if (from_file == true){
      return ERROR_ILLEGAL_ARGUMENT_COMBINATION;
}else{
      show_coordinates = false;
}
    }

  if (!embedded_format.empty() && !(output_format == "sdf" && (embedded_format == "inchi" || embedded_format == "smi"
                                    || embedded_format == "can")))
    {
      std::cerr << "Embedded format option is only possible if output format is SDF and option can have only inchi, smi, or can values." << std::endl;
      return ERROR_ILLEGAL_ARGUMENT_COMBINATION;
    }

  // This will hide the output "Warning: non-positive median line gap" from GOCR. Remove after this is fixed:
  fclose(stderr);
  OpenBabel::obErrorLog.StopLogging();

  bool is_reaction = false;
  if (output_format == "cmlr" || output_format == "rsmi" || output_format =="rxn")
    is_reaction = true;


  std::vector<std::vector<std::string> > pages_of_structures(page, std::vector<std::string> (0));
  std::vector<std::vector<Image> > pages_of_images(page, std::vector<Image> (0));
  std::vector<std::vector<double> > pages_of_avg_bonds(page, std::vector<double> (0));
  std::vector<std::vector<double> > pages_of_ind_conf(page, std::vector<double> (0));
  std::vector<std::vector<box_t> > pages_of_boxes(page, std::vector<box_t> (0));
  std::vector<std::vector<arrow_t> > arrows(page, std::vector<arrow_t>(0));
  std::vector<std::vector<plus_t> > pluses(page, std::vector<plus_t>(0));

  int total_structure_count = 0;
  int num_resolutions = NUM_RESOLUTIONS;
  if (input_resolution != 0)
    num_resolutions = 1;
  std::vector<double> array_of_confidence(num_resolutions, 0);
  std::vector<int> boxes_per_res(num_resolutions,0);
  std::vector<int> select_resolution(num_resolutions, input_resolution);
  set_select_resolution(select_resolution,input_resolution);
  std::vector<std::vector<std::vector<std::string> > > array_of_structures_page(
      page, std::vector<std::vector<std::string> >(num_resolutions));
  std::vector<std::vector<std::vector<double> > > array_of_avg_bonds_page(
      page, std::vector<std::vector<double> >(num_resolutions));
  std::vector<std::vector<std::vector<double> > > array_of_ind_conf_page(
      page, std::vector<std::vector<double> >(num_resolutions));
  std::vector<std::vector<std::vector<Image> > > array_of_images_page(
      page, std::vector<std::vector<Image> > (num_resolutions));
  std::vector<std::vector<std::vector<box_t> > > array_of_boxes_page(
      page, std::vector<std::vector<box_t> >(num_resolutions));

//#pragma omp parallel for default(shared) private(OCR_JOB,JOB)
  for (int l = 0; l < page; l++)
    {
      Image image;
      double page_scale=1;
      poppler::page_renderer poppler_renderer;

      int ttt = 0;

      if (verbose)
        std::cout << "Processing page " << (l+1) << " out of " << page << "..." << std::endl;

      if (type == "PDF" || type == "PS")
        page_scale *= (double) 72 / input_resolution;


      if (poppler_doc) // process PDF and PS files
	{
	  int resolution = input_resolution;
	  if (resolution == 0)
	    resolution = 300;
	  image = process_pdf_page(poppler_doc, poppler_renderer, l, resolution);
	}
      else
	{

          std::ostringstream pname;
	  pname << input_file << "[" << l << "]";
#pragma omp critical
	  {
	    image.read(pname.str());
	  }

	}
      if (l == 0 && !preview.empty())
	{
	  try
	    {
	      image.write(preview);
	    }
	  catch(...)
	    {}
	}

      image.modifyImage();
      bool adaptive = convert_to_gray(image, invert, adaptive_option, verbose);

      std::vector<std::vector<std::string> > array_of_structures(num_resolutions);
      std::vector<std::vector<double> > array_of_avg_bonds(num_resolutions), array_of_ind_conf(num_resolutions);
      std::vector<std::vector<Image> > array_of_images(num_resolutions);
      std::vector<std::vector<box_t> > array_of_boxes(num_resolutions);

      if (input_resolution > 300)
        {
          int percent = (100 * 300) / input_resolution;
          std::ostringstream scale;
          scale << percent << "%";
          image.scale(scale.str());
          page_scale /= (double) percent / 100;
        }

      if (verbose)
        std::cout << "Input resolutions are " << select_resolution << std::endl;

      ColorGray bgColor = getBgColor(image);
      if (rotate != 0)
        {
          image.backgroundColor(bgColor);
          image.rotate(rotate);
        }

      double rotation = 0;
      int unpaper_dx = 0;
      int unpaper_dy = 0;
      for (int i = 0; i < do_unpaper; i++)
        {
          double radians=0;
          int dx=0, dy=0;
          unpaper(image,radians,dx,dy);
          rotation +=radians;
          unpaper_dx +=dx;
          unpaper_dy +=dy;
        }

      // 0.1 is used for THRESHOLD_BOND here to allow for farther processing.
      std::list<std::list<std::list<point_t> > > clusters = find_segments(image, 0.1, bgColor, adaptive, is_reaction, arrows[l], pluses[l], verbose);

      if (verbose)
        std::cout << "Number of clusters: " << clusters.size() << '.' << std::endl;

      std::vector<box_t> boxes;
      std::set<std::pair<int, int> > brackets;
      int n_boxes = prune_clusters(clusters, boxes, brackets);
      std::sort(boxes.begin(), boxes.end(), comp_boxes);

      if (verbose)
        std::cout << "Number of boxes: " << boxes.size() << '.' << std::endl;

      // Creating r-group maps for all possible combos

        std::vector<std::map<std::string, std::string> > list_of_rgroup_maps;

        // Instatiate and populate map and vector
        std::map<std::string, std::string> map_of_rgroups;
        std::map<std::string, std::string> dummy_map_of_rgroups;
        std::vector<std::string> rgroup_vars;

        map_of_rgroups.insert(std::make_pair("R", "CH3"));
        dummy_map_of_rgroups.insert(std::make_pair("R", "OCH3"));

        rgroup_vars.push_back("R");
        list_of_rgroup_maps.push_back(map_of_rgroups);
        list_of_rgroup_maps.push_back(dummy_map_of_rgroups);

        for (int q = 0; q < list_of_rgroup_maps.size(); q++) {

            for (int res_iter = 0; res_iter < num_resolutions; res_iter++) {
                int total_boxes = 0;
                double total_confidence = 0;

                int resolution = select_resolution[res_iter];
                int working_resolution = resolution;
                if (resolution > 300)
                    working_resolution = 300;

                double THRESHOLD_BOND = set_threshold(threshold, resolution);

                int max_font_height = MAX_FONT_HEIGHT * working_resolution / 150;
                int max_font_width = MAX_FONT_WIDTH * working_resolution / 150;
                bool thick = true;
                if (resolution < 150)
                    thick = false;
                else if (resolution == 150 && !jaggy)
                    thick = false;

                //Image dbg = image;
                //dbg.modifyImage();
                //dbg.backgroundColor("white");
                //dbg.erase();
                //dbg.type(TrueColorType);
                for (int k = 0; k < n_boxes; k++)
                    if ((boxes[k].x2 - boxes[k].x1) > max_font_width && (boxes[k].y2 - boxes[k].y1) > max_font_height
                        && !boxes[k].c.empty() && ((boxes[k].x2 - boxes[k].x1) > 2 * max_font_width || (boxes[k].y2
                                                                                                        - boxes[k].y1) >
                                                                                                       2 *
                                                                                                       max_font_height)) {
                        int n_atom = 0, n_bond = 0, n_letters = 0, n_label = 0;
                        std::vector<atom_t> atom;
                        std::vector<bond_t> bond;
                        std::vector<letters_t> letters;
                        std::vector<label_t> label;
                        double box_scale = 1;
                        Image orig_box(Geometry(boxes[k].x2 - boxes[k].x1 + 2 * FRAME, boxes[k].y2 - boxes[k].y1 + 2
                                                                                                                   *
                                                                                                                   FRAME),
                                       bgColor);

                        for (unsigned int p = 0; p < boxes[k].c.size(); p++) {
                            int x = boxes[k].c[p].x;
                            int y = boxes[k].c[p].y;
                            ColorGray color = image.pixelColor(x, y);
                            //dbg.pixelColor(x, y, color);
                            orig_box.pixelColor(x - boxes[k].x1 + FRAME, y - boxes[k].y1 + FRAME, color);
                        }


                        int width = orig_box.columns();
                        int height = orig_box.rows();
                        Image thick_box;
                        create_thick_box(orig_box, thick_box, width, height, resolution, working_resolution, box_scale,
                                         bgColor, THRESHOLD_BOND, res_iter, thick, jaggy);

                        if (verbose)
                            std::cout << "Analysing box " << boxes[k].x1 << "x" << boxes[k].y1 << "-" << boxes[k].x2
                                      << "x" << boxes[k].y2 << " using working resolution " << working_resolution << '.'
                                      << std::endl;

                        Image box;
                        if (thick)
                            box = thin_image(thick_box, THRESHOLD_BOND, bgColor);
                        else
                            box = thick_box;
                        potrace_state_t *const st = raster_to_vector(box, bgColor, THRESHOLD_BOND, width, height,
                                                                     working_resolution);
                        potrace_path_t const *const p = st->plist;
                        n_atom = find_atoms(p, atom, bond, &n_bond, width, height);

                        int real_font_width, real_font_height;
                        n_letters = find_chars_rgroup(p, orig_box, letters, atom, bond, n_atom, n_bond, height, width,
                                                      bgColor,
                                                      THRESHOLD_BOND, max_font_width, max_font_height, real_font_width,
                                                      real_font_height, verbose, "RX");
                        if (verbose)
                            std::cout << "Number of atoms: " << n_atom << ", bonds: " << n_bond << ", " << n_letters
                                      << " letters: " << n_letters << " " << letters << " after find_atoms()"
                                      << std::endl;

                        double avg_bond_length = percentile75(bond, n_bond, atom);

                        double max_area = avg_bond_length * 5;
                        if (thick)
                            max_area = avg_bond_length;
                        n_letters = find_plus_minus(p, orig_box, bgColor, THRESHOLD_BOND, letters, atom, bond, n_atom,
                                                    n_bond, height, width,
                                                    real_font_height, real_font_width, n_letters, avg_bond_length);
                        n_atom = find_small_bonds(p, atom, bond, n_atom, &n_bond, max_area, avg_bond_length / 2, 5);

                        //remove_small_bonds_in_chars(atom,bond,letters);

                        find_old_aromatic_bonds(p, bond, n_bond, atom, n_atom, avg_bond_length);

                        if (verbose)
                            std::cout << "Number of atoms: " << n_atom << ", bonds: " << n_bond << ", " << n_letters
                                      << "letters: " << letters << " after find_old_aromatic_bonds()" << std::endl;

                        double dist = 3.;
                        if (working_resolution < 150)
                            dist = 2;

                        double thickness = skeletize(atom, bond, n_bond, box, THRESHOLD_BOND, bgColor, dist,
                                                     avg_bond_length);
                        remove_disconnected_atoms(atom, bond, n_atom, n_bond);
                        collapse_atoms(atom, bond, n_atom, n_bond, 3);
                        remove_zero_bonds(bond, n_bond, atom);

                        n_bond = find_wavy_bonds(bond, n_bond, atom, avg_bond_length);
                        //				if (ttt++ == 0)  debug_image(orig_box, atom, n_atom, bond, n_bond, "tmp.png");
                        n_letters = find_fused_chars(bond, n_bond, atom, letters, n_letters, real_font_height,
                                                     real_font_width, 0, orig_box, bgColor, THRESHOLD_BOND, 3, verbose);

                        n_letters = find_fused_chars(bond, n_bond, atom, letters, n_letters, real_font_height,
                                                     real_font_width, '*', orig_box, bgColor, THRESHOLD_BOND, 5,
                                                     verbose);

                        flatten_bonds(bond, n_bond, atom, 3);
                        remove_zero_bonds(bond, n_bond, atom);
                        avg_bond_length = percentile75(bond, n_bond, atom);

                        if (verbose)
                            std::cout << "Average bond length: " << avg_bond_length << std::endl;

                        double max_dist_double_bond = dist_double_bonds(atom, bond, n_bond, avg_bond_length);
                        n_bond = double_triple_bonds(atom, bond, n_bond, avg_bond_length, n_atom, max_dist_double_bond);
                        n_atom = find_dashed_bonds(p, atom, bond, n_atom, &n_bond,
                                                   std::max(MAX_DASH, int(avg_bond_length / 3)),
                                                   avg_bond_length, orig_box, bgColor, THRESHOLD_BOND, thick,
                                                   avg_bond_length, letters);

                        n_letters = remove_small_bonds(bond, n_bond, atom, letters, n_letters, real_font_height,
                                                       MIN_FONT_HEIGHT, avg_bond_length);

                        n_letters = find_numbers(p, orig_box, letters, atom, bond, n_atom, n_bond, height, width,
                                                 bgColor,
                                                 THRESHOLD_BOND, n_letters);

                        dist = 4.;
                        if (working_resolution < 300)
                            dist = 3;
                        if (working_resolution < 150)
                            dist = 2;

                        n_bond = fix_one_sided_bonds(bond, n_bond, atom, dist, avg_bond_length);

                        n_letters = clean_unrecognized_characters(bond, n_bond, atom, real_font_height, real_font_width,
                                                                  4,
                                                                  letters, n_letters);

                        thickness = find_wedge_bonds(thick_box, atom, n_atom, bond, n_bond, bgColor, THRESHOLD_BOND,
                                                     max_dist_double_bond, avg_bond_length, 3, 1);

                        n_label = assemble_labels(letters, n_letters, label);

                        if (verbose)
                            std::cout << n_label << " labels: " << label << " after assemble_labels()" << std::endl;

                        remove_disconnected_atoms(atom, bond, n_atom, n_bond);

                        collapse_atoms(atom, bond, n_atom, n_bond, thickness);

                        remove_zero_bonds(bond, n_bond, atom);

                        flatten_bonds(bond, n_bond, atom, 2 * thickness);

                        remove_zero_bonds(bond, n_bond, atom);

                        avg_bond_length = percentile75(bond, n_bond, atom);

                        collapse_double_bonds(bond, n_bond, atom, max_dist_double_bond);

                        extend_terminal_bond_to_label(atom, letters, n_letters, bond, n_bond, label, n_label,
                                                      avg_bond_length / 2,
                                                      thickness, max_dist_double_bond);

                        remove_disconnected_atoms(atom, bond, n_atom, n_bond);
                        collapse_atoms(atom, bond, n_atom, n_bond, thickness);
                        collapse_doubleup_bonds(bond, n_bond);

                        remove_zero_bonds(bond, n_bond, atom);
                        flatten_bonds(bond, n_bond, atom, thickness);
                        remove_zero_bonds(bond, n_bond, atom);
                        remove_disconnected_atoms(atom, bond, n_atom, n_bond);

                        extend_terminal_bond_to_bonds(atom, bond, n_bond, avg_bond_length, 2 * thickness,
                                                      max_dist_double_bond);

                        std::vector<bracket_t> bracket_boxes;
                        remove_bracket_atoms(atom, n_atom, bond, n_bond, brackets, thickness, boxes[k].x1, boxes[k].y1,
                                             box_scale, real_font_width, real_font_height, bracket_boxes);
                        remove_zero_bonds(bond, n_bond, atom);
                        remove_vertical_bonds_close_to_brackets(bracket_boxes, atom, bond, n_bond, thickness,
                                                                avg_bond_length);
                        remove_zero_bonds(bond, n_bond, atom);
                        flatten_bonds(bond, n_bond, atom, 2 * thickness);
                        assign_labels_to_brackets(bracket_boxes, label, n_label, letters, n_letters, real_font_width,
                                                  real_font_height);

                        collapse_atoms(atom, bond, n_atom, n_bond, 3);

                        remove_zero_bonds(bond, n_bond, atom);
                        flatten_bonds(bond, n_bond, atom, 5);
                        remove_zero_bonds(bond, n_bond, atom);

                        n_letters = clean_unrecognized_characters(bond, n_bond, atom, real_font_height, real_font_width,
                                                                  0,
                                                                  letters, n_letters);
                        int recognized_chars = count_recognized_chars(atom, bond);

                        std::map<std::string, std::string> current_rgroup = list_of_rgroup_maps[q];
                        std::cout << "Current R-Group R value : " << current_rgroup["R"] << std::endl;

                        for (int m = 0; m < n_atom; m++) {
//                      std::cout << "This atom is: " << atom[m] << std::endl;
//                      std::cout << "This atom's label is: " << atom[m].label << std::endl;
                            for (int z = 0; z < rgroup_vars.size(); z++) {


                                if (atom[m].label == rgroup_vars[z]) {
                                    std::cout << "Previous Atom was " << atom[m] << std::endl;

                                    atom[m].label = current_rgroup[rgroup_vars[z]];
                                    std::cout << " Atom updated to " << atom[m] << std::endl;
                                }

                            }

//                      if(atom[m].label == "R"){
//                          std::cout << "Previous Atom was " << atom[m] << std::endl;
//
//                          atom[m].label = "CH3";
//                          std::cout << " Atom updated to " << atom[m] << std::endl;
//                      }
                        }

                        assign_charge(atom, bond, n_atom, n_bond, spelling, superatom, debug);
                        find_up_down_bonds(bond, n_bond, atom, thickness);
                        int real_atoms = count_atoms(atom, n_atom);
                        int bond_max_type = 0;
                        int real_bonds = count_bonds(bond, n_bond, bond_max_type);

                        if (verbose)
                            std::cout << "Final number of atoms: " << real_atoms << ", bonds: " << real_bonds
                                      << ", chars: " << n_letters << '.' << std::endl;

                        std::cout << "Raw extracted atoms:  " << atom.size() << std::endl;
                        std::cout << "No of atoms:  " << n_atom << std::endl;
                        //std::cout << "All bonds:  " << bond << std::endl;


                        // std::cout << "No of bonds:  " << n_atom << std::endl;
//                  std::cout << "superatoms:  " << std::endl;
//                  for(std::map<string, pair<string,string> >::const_iterator it = superatom.begin();
//                          it != superatom.end(); ++it)
//                  {
//                      std::cout << it->first << " " << it->second.first << " " << it->second.second << "\n";
//                  }
                        //std::cout << "Real atoms:  " << n_atom << std::endl;


                        split_fragments_and_assemble_structure_record(atom, n_atom, bond, n_bond, boxes,
                                                                      l, k, resolution, res_iter,
                                                                      output_image_file_prefix, image, orig_box,
                                                                      real_font_width, real_font_height,
                                                                      thickness, avg_bond_length, superatom, real_atoms,
                                                                      real_bonds, bond_max_type,
                                                                      box_scale, page_scale, rotation, unpaper_dx,
                                                                      unpaper_dy, output_format, embedded_format,
                                                                      is_reaction, show_confidence,
                                                                      show_resolution_guess, show_page,
                                                                      show_coordinates,
                                                                      show_avg_bond_length, array_of_structures,
                                                                      array_of_avg_bonds, array_of_ind_conf,
                                                                      array_of_images, array_of_boxes, total_boxes,
                                                                      total_confidence,
                                                                      recognized_chars, show_learning, res_iter,
                                                                      verbose,
                                                                      bracket_boxes);

                        if (st != NULL)
                            potrace_state_free(st);

                    }
                array_of_confidence[res_iter] += total_confidence;
                boxes_per_res[res_iter] += total_boxes;
                //dbg.write("debug.png");
            }
        }

      #pragma omp critical
      {
         if (show_learning)
	  for (int j = 0; j < num_resolutions; j++)
	    for (unsigned int i = 0; i < array_of_structures[j].size(); i++)
	    {
	      pages_of_structures[l].push_back(array_of_structures[j][i]);
	      if (!output_image_file_prefix.empty())
		pages_of_images[l].push_back(array_of_images[j][i]);
	      pages_of_avg_bonds[l].push_back(array_of_avg_bonds[j][i]);
	      pages_of_ind_conf[l].push_back(array_of_ind_conf[j][i]);
	      pages_of_boxes[l].push_back(array_of_boxes[j][i]);
	      total_structure_count++;
	    }
	 else
	   for (int j = 0; j < num_resolutions; j++)
	     {
                array_of_structures_page[l][j] = array_of_structures[j];
		if (!output_image_file_prefix.empty())
		  array_of_images_page[l][j] = array_of_images[j];
		array_of_avg_bonds_page[l][j] = array_of_avg_bonds[j];
		array_of_ind_conf_page[l][j] = array_of_ind_conf[j];
		array_of_boxes_page[l][j] = array_of_boxes[j];
	      }

       }
     }

    double max_conf = -FLT_MAX;
    int max_res = 0;
    for (int i = 0; i < num_resolutions; i++)
        {
          if (boxes_per_res[i] > 0 && array_of_confidence[i]/boxes_per_res[i] > max_conf)
            {
              max_conf = array_of_confidence[i]/boxes_per_res[i];
              max_res = i;
            }
        }
      for (int i = 0; i < num_resolutions; i++)
	if (boxes_per_res[i] > 0 && array_of_confidence[i]/boxes_per_res[i] == max_conf && select_resolution[i] == 300) // second 300 dpi is without thinning
	  {
	    max_res = i;
	    break;
	  }

	if (!show_learning)
         for (int l = 0; l < page; l++)
	    {
	      pages_of_structures[l] = array_of_structures_page[l][max_res];
	      if (!output_image_file_prefix.empty())
		pages_of_images[l] = array_of_images_page[l][max_res];
	      pages_of_avg_bonds[l] = array_of_avg_bonds_page[l][max_res];
	      pages_of_ind_conf[l] = array_of_ind_conf_page[l][max_res];
	      pages_of_boxes[l] = array_of_boxes_page[l][max_res];
	      total_structure_count += array_of_structures_page[l][max_res].size();
	    }

  double best_bond = 0;

  //if (total_structure_count >= STRUCTURE_COUNT)
  //  find_limits_on_avg_bond(best_bond, pages_of_avg_bonds, pages_of_ind_conf);

  // If multiple pages are processed at several  resolutions different pages
  // may be processed at different resolutions leading to a seemingly different average bond length
  // Currently multi-page documents (PDF and PS) are all processed at the same resolution
  // and single-page images have all structures on the page at the same resolution

  //cout << min_bond << " " << max_bond << endl;


std::ostream &out_stream = outfile.is_open() ? outfile : std::cout;


  // For Andriod version we will find the structure with maximum confidence value, as the common usecase for Andriod is to analyse the
  // image (taken by embedded photo camera) that usually contains just one molecule:
  double max_confidence = -FLT_MAX;
  int l_index = 0;
  int i_index = 0;
  int image_count = 0;

  for (int l = 0; l < page; l++)
    {
      for (unsigned int i = 0; i < pages_of_structures[l].size(); i++)
	if (best_bond == 0 || (pages_of_avg_bonds[l][i] > best_bond/2 && pages_of_avg_bonds[l][i] < 2*best_bond))
	  {
	    if (pages_of_ind_conf[l][i] > max_confidence)
	      {
		max_confidence = pages_of_ind_conf[l][i];
		l_index = l;
		i_index = i;
	      }

	    if (output_format != "mol" && !is_reaction)
	      {
		out_stream << pages_of_structures[l][i];

		// Dump this structure into a separate file:
		if (!output_image_file_prefix.empty())
		  {
                    std::ostringstream fname;
		    fname << output_image_file_prefix << image_count << ".png";
		    image_count++;
		    if (fname.str() != "")
		      {
			Image tmp = pages_of_images[l][i];
			if (resize != "")
			  {
			    tmp.scale(resize);
			  }
			tmp.write(fname.str());
		      }
		  }
	      }
	  }
      if (is_reaction && !arrows[l].empty())
	 {
           std::vector<std::string> reactions;
           std::vector<box_t> rbox;
	   arrange_reactions(arrows[l], pages_of_boxes[l], pluses[l], reactions, rbox, pages_of_structures[l],output_format);
	   for (int k=0; k<reactions.size(); k++)
	     {
               out_stream << reactions[k] << std::endl;

	       if (!output_image_file_prefix.empty())
		 {
                   std::ostringstream fname;
		   fname << output_image_file_prefix << image_count << ".png";
		   image_count++;
		   if (fname.str() != "")
		     {
		       Image tmp = pages_of_images[l][k];
		       Geometry geometry = Geometry(rbox[k].x2 - rbox[k].x1, rbox[k].y2 - rbox[k].y1, rbox[k].x1, rbox[k].y1);
		       tmp.crop(geometry);
		       if (resize != "")
			 {
			   tmp.scale(resize);
			 }
		       tmp.write(fname.str());
		     }
		 }
	     }
	 }
    }
  // Output the structure with maximum confidence value:
  if (output_format == "mol")
    {
      out_stream << pages_of_structures[l_index][i_index];
      if (!output_image_file_prefix.empty())
	{
          std::ostringstream fname;
	  fname << output_image_file_prefix  << ".png";
	  if (fname.str() != "")
	    {
	      Image tmp = pages_of_images[l_index][i_index];
	      if (resize != "")
		{
		  tmp.scale(resize);
		}
	      tmp.write(fname.str());
	    }
	}
    }

  out_stream.flush();

if (from_file != true){
  if (!output_file.empty())
    outfile.close();
}

  return 0;
}

void test_osra_lib(const std::string &output, int pointless)
{
  using namespace std;

  std::cout << output << std::endl;

#ifdef OSRA_LIB
  std::cout << "osra lib found" << endl;
#else
  std::cout << "osra lib NOT found" << endl;
#endif

  return;
}