
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>

#include "nvutility.h"

#include "bag.h"
#include "chrtr2.h"

#include "version.h"


/*

  bag_residual

  Jan C. Depner

  July 14th 2010

  This program compares two BAG files and creates a CHRTR2 difference surface.
  The files must have the same extents and bin sizes.

*/

#define SKOSH 0.00000001L

void usage ()
{
  fprintf (stderr, "\nUsage: bag_residual BAG1 BAG2\n");
  fprintf (stderr, "\n\n");
  fprintf (stderr, "\tBAG2 will be subtracted from BAG1.  A CHRTR2 file of the\n");
  fprintf (stderr, "\tdifference surface will be created.  The file will be\n");
  fprintf (stderr, "\tnamed BAG1.ch2.\n\n");
  exit (-1);
}


int32_t main (int32_t argc, char **argv)
{
  bagHandle       bagHandle[2];
  bagError        err;
  u8              *errstr;
  int32_t         i, j, percent, old_percent, neg_count, pos_count, option_index, chrtr2_handle = -1,
                  bin_width[2], bin_height[2], total_bins;
  float           *data[2], rms, max_val, min_val, neg_percent, pos_percent, depthtot, min_depth, max_depth, diff = 0.0;
  double          x_bin_size_degrees[2], y_bin_size_degrees[2], sum, sum2, meandiff, meandepth, ss, var, stddev, sddepth, dep = 0.0;
  NV_F64_XYMBR    mbr[2];
  CHRTR2_HEADER   chrtr2_header;
  CHRTR2_RECORD   chrtr2_record;
  char            bag_file[2][512], chrtr2_file[512];
  char            c;
  extern char     *optarg;
  extern int      optind;


  fprintf(stderr, "\n\n %s \n\n", VERSION);
  fflush (stderr);


  option_index = 0;
  while (NVTrue) 
    {
      static struct option long_options[] = {{0, no_argument, 0, 0}};

      c = (char) getopt_long (argc, argv, "b", long_options, &option_index);
      if (c == -1) break;

      switch (c) 
        {
        case 0:

          switch (option_index)
            {
            case 0:
              break;
            }
          break;

        case 'b':
          /*  Place holder  */
          break;

        default:
          usage ();
          break;
        }
    }


  /* Make sure we got the mandatory file names.  */
  
  if (optind >= argc) usage ();


  /*  Override the HDF5 version check so that we can read BAGs created with an older version of HDF5.  */

  putenv ("HDF5_DISABLE_VERSION_CHECK=2");


  strcpy (bag_file[0], argv[optind]);
  strcpy (bag_file[1], argv[optind + 1]);


  for (i = 0 ; i < 2 ; i++)
    {
      if ((err = bagFileOpen (&bagHandle[i], BAG_OPEN_READONLY, (u8 *) bag_file[i])) != BAG_SUCCESS)
        {
          fprintf (stderr, "\nError opening BAG file %s\n", bag_file[i]);
          if (bagGetErrorString (err, &errstr) == BAG_SUCCESS) fprintf (stderr, "%s\n", errstr);
          exit (-1);
        }


      bin_width[i] = bagGetDataPointer (bagHandle[i])->def.ncols;
      bin_height[i] = bagGetDataPointer (bagHandle[i])->def.nrows;
      x_bin_size_degrees[i] = bagGetDataPointer (bagHandle[i])->def.nodeSpacingX;
      y_bin_size_degrees[i] = bagGetDataPointer (bagHandle[i])->def.nodeSpacingY;
      mbr[i].min_x = bagGetDataPointer (bagHandle[i])->def.swCornerX;
      mbr[i].min_y = bagGetDataPointer (bagHandle[i])->def.swCornerY;
      mbr[i].max_x = mbr[i].min_x + bin_width[i] * x_bin_size_degrees[i];
      mbr[i].max_y = mbr[i].min_y + bin_height[i] * y_bin_size_degrees[i];


      data[i] = (float *) calloc (bin_width[i], sizeof (float));

      if (data[i] == NULL)
        {
          perror ("Allocating data row");
          exit (-1);
        }
    }


  /*  Check for (near) identical extents and spacing.  */

  if (abs (bin_width[0] - bin_width[1]) > 2 || abs (bin_height[0] - bin_height[1]) > 2 || fabs (mbr[0].min_x - mbr[1].min_x) > SKOSH ||
      fabs (mbr[0].min_y - mbr[1].min_y) > SKOSH || fabs (x_bin_size_degrees[0] - x_bin_size_degrees[1]) > SKOSH ||
      fabs (y_bin_size_degrees[0] - y_bin_size_degrees[1]) > SKOSH)
    {
      fprintf (stderr, "\nBAG file extents and/or spacing do not match.\n");
      fprintf (stderr, "BAG1 MinX = %.7f  MinY = %.7f  X = %.7f  Y = %.7f\n", mbr[0].min_x, mbr[0].min_y, x_bin_size_degrees[0], y_bin_size_degrees[0]);
      fprintf (stderr, "BAG1 MinX = %.7f  MinY = %.7f  X = %.7f  Y = %.7f\n\n", mbr[1].min_x, mbr[1].min_y, x_bin_size_degrees[1], y_bin_size_degrees[1]);
      free (data[0]);
      free (data[1]);
      exit (-1);
    }


  memset (&chrtr2_header, 0, sizeof (CHRTR2_HEADER));


  /*  Generate the chrtr2 file name.  */

  strcpy (chrtr2_file, bag_file[0]);
  strcpy (&chrtr2_file[strlen (chrtr2_file) - 4], ".ch2");


  /*  Populate the chrtr2 header prior to creating the file.  */

  strcpy (chrtr2_header.creation_software, VERSION);
  chrtr2_header.z_units = CHRTR2_METERS;
  chrtr2_header.mbr.wlon = mbr[0].min_x;
  chrtr2_header.mbr.slat = mbr[0].min_y;
  chrtr2_header.width = bin_width[0];
  chrtr2_header.height = bin_height[0];
  chrtr2_header.lat_grid_size_degrees = y_bin_size_degrees[0];
  chrtr2_header.lon_grid_size_degrees = x_bin_size_degrees[0];
  chrtr2_header.min_z = -326.00;
  chrtr2_header.max_z = 326.00;
  chrtr2_header.z_scale = 100.0;
  chrtr2_header.horizontal_uncertainty_scale = 0.0;
  chrtr2_header.vertical_uncertainty_scale = 0.0;


  /*  Try to create and open the chrtr2 file.  */

  chrtr2_handle = chrtr2_create_file (chrtr2_file, &chrtr2_header);
  if (chrtr2_handle < 0)
    {
      chrtr2_perror ();
      exit (-1);
    }


  fprintf(stderr, "\n\n");
  fflush (stderr);


  percent = 0;
  old_percent = -1;


  total_bins = bin_height[0] * bin_width[0];


  neg_count = 0;
  pos_count = 0;
  min_val = 99999.0;
  max_val = -99999.0;
  min_depth = 99999.0;
  max_depth = -99999.0;
  sum = 0.0;
  sum2 = 0.0;
  depthtot = 0.0;


  /* Process all records in the BAG files */

  for (i = 0 ; i < bin_height[0] ; i++)
    {
      err = bagReadRow (bagHandle[0], i, 0, bin_width[0] - 1, Elevation, (void *) data[0]);
      err = bagReadRow (bagHandle[1], i, 0, bin_width[1] - 1, Elevation, (void *) data[1]);

      for (j = 0 ; j < bin_width[0] ; j++)
        {
	  if (data[0][j] < NULL_ELEVATION && data[1][j] < NULL_ELEVATION)
            {
              dep = data[0][j];
              diff = data[0][j] - data[1][j];


              memset (&chrtr2_record, 0, sizeof (CHRTR2_RECORD));

              chrtr2_record.z = diff;
              chrtr2_record.status = CHRTR2_REAL;

              chrtr2_write_record_row_col (chrtr2_handle, i, j, chrtr2_record);


              if (dep < min_depth) min_depth = dep;

              if (dep > max_depth) max_depth = dep;

              depthtot += dep;
              sum += diff;
              sum2 += diff * diff;

              if (diff < 0.0)
                {
                  neg_count++;
                }
              else
                {
                  pos_count++;
                }

              if (fabs((double) diff) < min_val) min_val = fabs((double) diff);

              if (fabs((double) diff) > max_val) max_val = fabs((double) diff);
            }

          percent = ((float) (i * bin_width[0] + j) / (float) total_bins) * 100.0;

          if (old_percent != percent)
            {
              fprintf (stderr, "%03d%% processed     \r", percent);
              fflush (stderr);
              old_percent = percent;
            }
        }
    }

  percent = 100;
  fprintf (stderr, "%03d%% processed        \n\n", percent);
  fflush (stderr);


  printf ("#FIRST BAG file  : %s\n", bag_file[0]);
  printf ("#SECOND BAG file : %s\n#\n", bag_file[1]);


  chrtr2_close_file (chrtr2_handle);


  printf ("#       RMS       MEAN DIFF          STD             STD%%    NEG%%   POS%%      MAX RESID    MEAN DEPTH    # POINTS\n#\n");

  meandiff = sum / (float) (neg_count + pos_count);
  meandepth = depthtot / (float) (neg_count + pos_count);
  ss = sum2 - (sum * meandiff);
  var = ss / ((neg_count + pos_count) - 1);
  stddev = sqrt (var);
  sddepth = (stddev / meandepth) * 100;
  rms = sqrt((double) (sum2 / (float) (neg_count + pos_count)));
  neg_percent = ((float) neg_count / (float) (neg_count + pos_count)) * 100.0;
  pos_percent = ((float) pos_count / (float) (neg_count + pos_count)) * 100.0;

  if (sum != 0.0)
    {
      printf (" %10.3f   %10.3f      %10.3f      %10.4f    %03d    %03d   %10.3f    %10.3f  %12d\n", 
              rms, meandiff, stddev, sddepth, NINT (neg_percent), NINT (pos_percent), 
              max_val, meandepth, (neg_count + pos_count));
    }

  printf ("\n\n\n");


  return (0);
}
