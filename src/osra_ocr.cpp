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

#include <string.h> // strlen(), memset()
#include <ctype.h> // isalnum(), isspace()

#include <vector> // std:vector
#include <iostream> // std::cout

#include "osra_common.h"
#include "config.h"

extern "C" {
#include <pgm2asc.h>
}

#include <ocradlib.h>

#include "osra.h"
#include "osra_ocr.h"

#ifdef HAVE_CUNEIFORM_LIB
#include <cuneiform.h>
#endif

#ifdef HAVE_TESSERACT_LIB
// We can't push these functions into this cpp file, as the types from Tesseract conflict with GOCR
void osra_tesseract_init();
void osra_tesseract_destroy();
char osra_tesseract_ocr(unsigned char *pixel_map, int width, int height, const std::string &char_filter);
#endif

// Global GOCR variable (omg) both for 0.48-0.49 and 0.50 versions:
job_t *JOB;
job_t *OCR_JOB;

// Also declared in osra_ocr_tesseract.cpp:
const char UNKNOWN_CHAR = '_';

/**
 * THRESHOLD is the graylevel binarization threshold.
 * dropx and dropy are the coordinates for the starting point from where the connected component (the image of the character) will be searched for.
 * Very often there is no bounding rectangle that would exclude all extraneous pieces without cutting into the character itself.
 * Those pieces confuse the OCR libraries quite a bit, so it's better to extract connected components (all characters that OSRA needs to resolve
 * are luckily single connected components) and leave the extra bits out.
 */

void osra_ocr_init()
{
#ifdef HAVE_CUNEIFORM_LIB
  // Initialization for Cuneiform library should be called only once. Otherwise it breaks down during deinitialization:
  int langcode = PUMA_LANG_ENGLISH;
  PumaBool dotmatrix = 0;
  PumaBool fax = 0;
  PumaBool onecolumn = 1;
  PUMA_Init(0, 0);
  PUMA_SetImportData(PUMA_Word32_Language, &langcode);
  PUMA_SetImportData(PUMA_PumaBool32_DotMatrix, &dotmatrix);
  PUMA_SetImportData(PUMA_PumaBool32_Fax100, &fax);
  PUMA_SetImportData(PUMA_PumaBool32_OneColumn, &onecolumn);
#endif

#ifdef HAVE_TESSERACT_LIB
  osra_tesseract_init();
#endif
}

void osra_ocr_destroy()
{
#ifdef HAVE_CUNEIFORM_LIB
  PUMA_Done();
#endif

#ifdef HAVE_TESSERACT_LIB
  osra_tesseract_destroy();
#endif
}

// Function: osra_gocr_ocr()
//      Make an attempt to OCR the image box with GOCR engine.
//
// Parameters:
//      job_t - includes pixel map and character filter
//
// Returns:
//      0 in case the recognition failed or valid alphanumeric character
char osra_gocr_ocr(job_t &gocr_job)
{
  OCR_JOB = &gocr_job;
  JOB = &gocr_job;
// cout<<"Before gocr"<<endl;
  try
    {
      pgm2asc(&gocr_job);
    }
  catch (...)
    {
    }
//cout<<"After gocr"<<endl;
  char *l = (char *) gocr_job.res.linelist.start.next->data;

  if (l != NULL && strlen(l) == 1 && isalnum(l[0]))
    return l[0];

  return UNKNOWN_CHAR;
}

// Function: osra_ocrad_ocr()
//      Make an attempt to OCR the image box with OCRAD engine.
//
// Parameters:
//      ocrad_pixmap - includes pixel map and the image mode
//      char_filter - character filter
//
// Returns:
//      0 in case the recognition failed or valid alphanumeric character
char osra_ocrad_ocr(const OCRAD_Pixmap * const ocrad_pixmap, const std::string &char_filter)
{
  char result = 0;
  std::string line;

  OCRAD_Descriptor * const ocrad_res = OCRAD_open();

  // If the box height is less than 10px, it should be scaled up a bit, otherwise OCRAD is unable to catch it:
  if (ocrad_res && OCRAD_get_errno(ocrad_res) == OCRAD_ok && OCRAD_set_image(ocrad_res, ocrad_pixmap, 0) == 0
      && (ocrad_pixmap->height >= 10 || OCRAD_scale(ocrad_res, 2) == 0) && OCRAD_recognize(ocrad_res, 0) == 0)
    {
      result = OCRAD_result_first_character(ocrad_res);
      if (OCRAD_result_blocks(ocrad_res) >= 1 && OCRAD_result_lines(ocrad_res, 0) && OCRAD_result_line(
            ocrad_res, 0, 0) != 0)
        line = OCRAD_result_line(ocrad_res, 0, 0);
    }

  OCRAD_close(ocrad_res);

  // TODO: Why line should have 0 or 1 characters? Give examples...
  if (line.length() > 2 || !isalnum(result) || (!char_filter.empty() && char_filter.find(result, 0) == std::string::npos))
    return UNKNOWN_CHAR;

  return result;
}

// Function: osra_cuneiform_ocr()
//      Make an attempt to OCR the image box with Cuneiform engine.
//
// Parameters:
//      cuneiform_img - pixel map
//      verbose - if set, then output intermediate results
//      char_filter - character filter
//
// Returns:
//      0 in case the recognition failed or valid alphanumeric character
#ifdef HAVE_CUNEIFORM_LIB
char osra_cuneiform_ocr(Magick::Image &cuneiform_img, const std::string &char_filter)
{
  Magick::Blob blob;
  cuneiform_img.write(&blob, "DIB");
  size_t data_size = blob.length();
  char *dib = new char[data_size];
  memcpy(dib, blob.data(), data_size);

  char str[256];
  memset(str, 0, sizeof(str));

  if (!PUMA_XOpen(dib, NULL) || !PUMA_XFinalRecognition() || !PUMA_SaveToMemory(NULL, PUMA_TOTEXT, PUMA_CODE_ASCII, str, sizeof(str) - 1))
    {
      //if (verbose)
      //  cout << "Cuneiform recognition failed." << endl;

      PUMA_XClose();
      delete []dib;

      return UNKNOWN_CHAR;
    }

  PUMA_XClose();
  delete []dib;

  // As we have initialized the image with two identical samples, it is expected that they go in the string
  // one after another, or separated by space (e.g. "ZZ\n" or "Z Z\n").
  if (((str[0] == str[1] && isspace(str[2])) || (str[0] == str[2] && str[1] == ' ')) && isalnum(str[0])
      && (char_filter.empty() || char_filter.find(str[0], 0) != std::string::npos))
    return str[0];

  return UNKNOWN_CHAR;
}
#endif

char get_atom_label(const Magick::Image &image, const Magick::ColorGray &bg, int x1, int y1, int x2, int y2,
                    double THRESHOLD, int dropx, int dropy, bool no_filtering, bool verbose, bool numbers)
{
  char c = UNKNOWN_CHAR;

  const int width = x2 - x1 + 1;
  const int height = y2 - y1 + 1;

  unsigned char *pixmap = (unsigned char *) malloc(width * height);

  for (int i = y1; i <= y2; i++)
    for (int j = x1; j <= x2; j++)
      pixmap[(i - y1) * width + j - x1] = (unsigned char) (255 - 255 * get_pixel(image, bg, j, i, THRESHOLD));

  // Here we drop down from the top of the box, middle of x coordinate and extract connected component
  int t = 1;
  int y = dropy - y1 + 1;
  int x = dropx - x1;

  while ((t != 0) && (y < height))
    {
      t = pixmap[y * width + x];
      y++;
    }

  if (t != 0)
    {
      free(pixmap);
      return 0;
    }

  #pragma omp critical
  {
    y--;

    pixmap[y * width + x] = 2;

    std::list<int> cx;
    std::list<int> cy;

    cx.push_back(x);
    cy.push_back(y);

    while (!cx.empty())
      {
        x = cx.front();
        y = cy.front();
        cx.pop_front();
        cy.pop_front();
        pixmap[y * width + x] = 1;

        // this goes around 3x3 square touching the chosen pixel
        for (int i = x - 1; i < x + 2; i++)
          for (int j = y - 1; j < y + 2; j++)
            if (i < width && j < height && i >= 0 && j >= 0 && pixmap[j * width + i] == 0)
              {
                cx.push_back(i);
                cy.push_back(j);
                pixmap[j * width + i] = 2;
              }
      }

    // Flatten the bitmap. Note: the bitmap is inverted after this cycle (255 means "empty", 0 means "pixel").
    for (int i = 0; i < height; i++)
      for (int j = 0; j < width; j++)
        pixmap[i * width + j] = (pixmap[i * width + j] == 1 ? 0 : 255);

    job_t gocr_job;

    // The list of all characters, that can be recognised as atom label:
    std::string char_filter = RECOGNIZED_CHARS;
    if (numbers) char_filter = "1";
    if (no_filtering) char_filter.clear();

    job_init(&gocr_job);
    job_init_image(&gocr_job);

    //gocr_job.cfg.cs = 160;
    //gocr_job.cfg.certainty = 80;
    //gocr_job.cfg.dust_size = 1;
    gocr_job.src.p.x = width;
    gocr_job.src.p.y = height;
    gocr_job.src.p.bpp = 1;
    gocr_job.src.p.p = pixmap;
    if (char_filter.empty())
      gocr_job.cfg.cfilter = (char*)NULL;
    else
      gocr_job.cfg.cfilter = (char*) char_filter.c_str();

    struct OCRAD_Pixmap *ocrad_pixmap = new OCRAD_Pixmap();
    unsigned char *ocrad_bitmap = (unsigned char *) malloc(width * height);

    memset(ocrad_bitmap, 0, width * height);

    ocrad_pixmap->height = height;
    ocrad_pixmap->width = width;
    ocrad_pixmap->mode = OCRAD_bitmap;
    ocrad_pixmap->data = ocrad_bitmap;

    // Number of non-zero pixels on the bitmap, excluding the 1px border:
    int pixmap_pixels_count = 0;
    // Number of zero pixels on the bitmap, excluding the 1px border:
    int pixmap_zeros_count = 0;

    // The code below initialises "opix->data" buffer ("bitmap_data") for OCRAD from "tmp" buffer:
#ifdef HAVE_CUNEIFORM_LIB
    Magick::Image cuneiform_img(Magick::Geometry(2 * width + 2, height), "white");
    // From cuneiform_src/cli/cuneiform-cli.cpp::preprocess_image(Magick::Image&):168
    cuneiform_img.monochrome();
    cuneiform_img.type(Magick::BilevelType);
#endif
    for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
          {
            if (pixmap[y * width + x] == 0)
              {
                ocrad_bitmap[y * width + x] = 1;
#ifdef HAVE_CUNEIFORM_LIB
                // Draw two identical samples that follow one another. We do so because Cuneiform has difficulties in recognizing single characters:
                cuneiform_img.pixelColor(x, y, "black");
                cuneiform_img.pixelColor(x + width + 2, y, "black");
#endif
                if (x > 0 && x < width - 1 && y > 0 && y < height - 1)
                  pixmap_pixels_count++;
              }
            else if (x > 0 && x < width - 1 && y > 0 && y < height - 1)
              pixmap_zeros_count++;
          }
      }

    if (verbose)
      {
        std::cout << "Box to OCR: " << x1 << "x" << y1 << "-" << x2 << "x" << y2 << " w/h: " << width << "x" << height << std::endl;
        for (int i = 0; i < height; i++)
          {
            for (int j = 0; j < width; j++)
              std::cout << (gocr_job.src.p.p[i * width + j] / 255 ? '#' : '.');
            std::cout << std::endl;
          }
      }

    if (pixmap_pixels_count <= MIN_CHAR_POINTS || pixmap_zeros_count <= MIN_CHAR_POINTS)
      goto FINALIZE;

    c = osra_gocr_ocr(gocr_job);

    if (verbose)
      std::cout << "GOCR: c=" << c << std::endl;

    //c = UNKNOWN_CHAR; // Switch off GOCR recognition

    // Character recognition succeeded for GOCR:
    if (c != UNKNOWN_CHAR)
      goto FINALIZE;

    // Character recognition failed for GOCR and we try OCRAD:
    c = osra_ocrad_ocr(ocrad_pixmap, char_filter);

    if (verbose)
      std::cout << "OCRAD: c=" << c << std::endl;

    //c = UNKNOWN_CHAR;  // Switch off OCRAD recognition

    // Character recognition succeeded for OCRAD:
    if (c != UNKNOWN_CHAR)
      goto FINALIZE;

#ifdef HAVE_TESSERACT_LIB
    c = osra_tesseract_ocr(gocr_job.src.p.p, width, height, char_filter);

    if (verbose)
      std::cout << "Tesseract: c=" << c << std::endl;

    //c = UNKNOWN_CHAR;  // Switch off Tesseract recognition

    // Character recognition succeeded for Tesseract:
    if (c != UNKNOWN_CHAR)
      goto FINALIZE;
#endif
#ifdef HAVE_CUNEIFORM_LIB
    // TODO: Why box width should be more than 7 for Cuneiform?
    if (width <= 7)
      goto FINALIZE;

    c = osra_cuneiform_ocr(cuneiform_img, char_filter);

    if (verbose)
      std::cout << "Cuneiform: c=" << c << std::endl;

    //c = UNKNOWN_CHAR; // Switch off Cuneiform recognition
#endif

FINALIZE:
    // "pixmap" is freed together with "gocr_job".
    job_free_image(&gocr_job);
    OCR_JOB = NULL;
    JOB = NULL;

    delete ocrad_pixmap; // delete OCRAD Pixmap
    free(ocrad_bitmap);

    // TODO: Why there are problems with "7" with a given box size? If the problem is engine-specific, it should be moved to appropriate section
    if (c == '7' && (width <= 10 || height <= 20))
      c = UNKNOWN_CHAR;
  } // #pragma omp critical

  return(c == UNKNOWN_CHAR ? 0 : c);
}

bool detect_square_bracket(unsigned char *pic, int x, int y)
{

  int w = -1;
  for (int j = x / 2; j >= 0 ; j--)
    {
      if (pic[(y / 2) * x + j] == 1)
	{
	  w = j;
	  break;
	}
    }
  int h1 = -1;
  for (int i = y / 2; i >= 0; i--)
    {
      if (pic[i * x + x / 2] == 1)
	{
	  h1 = i;
	  break;
	}
    }
  int h2 = y;
  for (int i = y / 2; i < y; i++)
    {
      if (pic[i * x + x / 2] == 1)
	{
	  h2 = i;
	  break;
	}
    }
  if (w > x /2 + 1 || h1 > y / 4 || (y - h2) > y / 4)
    return false;
  int total_vert(0), total_up_hor(0), total_down_hor(0), fill_vert(0), fill_up_hor(0), fill_down_hor(0), total_empty(0), fill_empty(0);

  for (int i = 0; i < y; i++)
    {
      for (int j = 0; j < x; j++)
	{
	  bool fill = (pic[i * x + j] == 1);
	  if (j <= w)
	    {
	      total_vert++;
	      if (fill)
		fill_vert++;
	    }
	  if (i <= h1)
	    {
	      total_up_hor++;
	      if (fill)
		fill_up_hor++;
	    }
	  if (i >= h2)
	    {
	      total_down_hor++;
	      if (fill)
		fill_down_hor++;
	    }
	  if (j > w && i > h1 && i < h2)
	    {
	      total_empty++;
	      if (fill)
		fill_empty++;
	    }
	}
    }

  if (total_vert == 0 || total_up_hor == 0 || total_down_hor == 0 || total_empty == 0)
    return false;

  double fill_v = double(fill_vert) / total_vert;
  double fill_up_h = double(fill_up_hor) / total_up_hor;
  double fill_down_h = double(fill_down_hor) / total_down_hor;
  double fill_e = double(fill_empty) / total_empty;
  if (fill_v > 0.8 && fill_up_h > 0.7 && fill_down_h > 0.7 && fill_e < 0.2)
    return true;
  return false;
}

bool detect_bracket(int x, int y, unsigned char *pic)
{
  bool res = false;


#pragma omp critical
  {
  char c1 = 0;
  job_t job;
  JOB = &job;
  OCR_JOB = &job;
  job_init(&job);
  job_init_image(&job);
  job.cfg.cfilter = (char *) "([{";

  //job.cfg.cs = 160;
  //job.cfg.certainty = 80;
  //job.cfg.dust_size = 1;

  job.src.p.x = x;
  job.src.p.y = y;
  job.src.p.bpp = 1;
  job.src.p.p = pic;

  struct OCRAD_Pixmap *ocrad_pixmap = new OCRAD_Pixmap();
  unsigned char *ocrad_bitmap = (unsigned char *) malloc(x * y);
  memset(ocrad_bitmap, 0, x * y);
  ocrad_pixmap->height = y;
  ocrad_pixmap->width = x;
  ocrad_pixmap->mode = OCRAD_bitmap;
  ocrad_pixmap->data = ocrad_bitmap;

  int count = 0;
  int zeros = 0;
  for (int i = 0; i < y; i++)
    {
     for (int j = 0; j < x; j++)
      {
	if (pic[i * x + j] == 0)
	  {
	    ocrad_bitmap[i * x + j] = 1;
	    count++;
	  }
	else
	  zeros++;
      }
   }

  /*std::string str;
  for (int i = 0; i < y; i++)
    {
      for (int j = 0; j < x; j++)
	{
	  str += (pic[i * x + j] != 0 ? "#" : ".");
	}
      str += "\n";
    }

    std::cout << str << std::endl;
  */
  if (count > MIN_CHAR_POINTS && zeros > MIN_CHAR_POINTS)
    {
      try
	{
	  pgm2asc(&job);
	}
      catch (...) {  }
      char *l;
      l = (char *) job.res.linelist.start.next->data;
      if (l != NULL)
	c1 = l[0];
      if (c1 == '(' || c1 == '[' || c1 == '{')
	res = true;
      else
	{
          char c2 = osra_ocrad_ocr(ocrad_pixmap, "([{");
	  if (c2 == '(' || c2 == '[' || c2 == '{')
	    res = true;
	}
    }

  if (!res)
    res = detect_square_bracket(ocrad_bitmap, x, y);
  /*
  if (res)
    {
      std::cout << "Found! " << c1 << std::endl;
    }
  */

  delete ocrad_pixmap;
  free(ocrad_bitmap);

  job_free_image(&job);
  OCR_JOB = NULL;
  JOB = NULL;
}

  return (res);
}


const std::string fix_atom_name(const std::string &s, int n,
                                const std::map<std::string, std::string> &fix,
                                const std::map<std::string, std::string> &superatom, bool debug)
{
  std::string r = s;

  if (s.length() == 1)
    r = toupper(s.at(0));
  if (s == "H" && n > 1)
    r = "N";

  std::map<std::string, std::string>::const_iterator it = fix.find(s);
  std::string mapped = " ";
  if (it != fix.end())
    {
      r = it->second;
      mapped = r;
    }

  if (debug && s != " " && s != "")
    {
      it = superatom.find(r);
      std::string smiles = " ";
      if (it != superatom.end())
        smiles = it->second;
      std::cout << s << " --> " << mapped << " --> " << smiles << std::endl;
    }

  return (r);
}
