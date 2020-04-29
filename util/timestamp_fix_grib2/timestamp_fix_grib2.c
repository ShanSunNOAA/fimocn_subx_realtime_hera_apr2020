#include <stdio.h>
#include <errno.h>
#include <regex.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <libgen.h>
#include <string.h>
#include <ctype.h>
#include "gribmsgfixtime.h"
#include "errmacs.h"

static int all_digits (char *s);

static const char USAGE[] =
  "The first argument must be a file name.  The file name may include\n"
  "the path to the file.  The second argument must be either of the form\n"
  "FFF (3 digits) or yyjjjHHMM0FFF (13 digits).  Calling this program\n"
  "should be like one of the following\n"
  "\n"
  "    timestamp_fix_grib2 <grib2 file>\n"
  "\n"
  "    timestamp_fix_grib2 <grib2 file> <FFF>\n"
  "\n"
  "    timestamp_fix_grib2 <grib2 file> <yyjjjHHMM0FFF>\n";

int
main (int argc, char *argv[])
{
  /* fn_hr meant file name hour but now kind o' means forecast time */
  uint64_t fn_hr = 0L;
  char *fname = NULL;

  /* get the value of fn_hr from the command line arguments */
  /* BEGIN{process the command line and arguments} */
  if (argc >= 2)
    {
      fname = argv[1];
    }

  switch (argc)
    {
    case 2:
      {
	char *bname;
	char *basec;
	char fnamepat[] = "\\([[:digit:]]\\{13\\}\\)";
	regex_t reg;
	int regerr;
	regmatch_t pmatch;

	errno = 0;

	regerr = regcomp (&reg, fnamepat, 0);
	if (regerr != 0)
	  {
	    char errbuff[80];
	    regerror (regerr, &reg, errbuff, sizeof errbuff);
	    fprintf (stderr, "%s\n", errbuff);
	    exit (1);
	  }
	ERREXIT ((basec = strdup (fname)) == NULL);
	bname = basename (basec);
	if (regexec (&reg, bname, 1, &pmatch, 0) == 0)
	  {
	    enum { timestamp_size = 14 }; /* 13 digits terminated by a '\0' */
	    char timestamp[timestamp_size];
	    ERREXIT (pmatch.rm_so < 0);
	    strncpy (timestamp, bname + pmatch.rm_so, timestamp_size);
	    if (ERRPRNT ((fn_hr = strtol (timestamp + 10, NULL, 10),
			  errno != 0)))
	      {
		exit (1);
	      }
	  }
	free (basec);
	regfree (&reg);
      }
      break;
    case 3:
      {
	size_t arg2_length = strlen (argv[2]);
	if (arg2_length == 3)
	  {
	    ERREXIT (! all_digits (argv[2]));
	    if (ERRPRNT ((fn_hr = (uint64_t) strtol (argv[2], NULL, 10),
			  errno != 0)))
	      {
		exit (1);
	      }
	  }
	else if (arg2_length == 13)
	  {
	    ERREXIT (! all_digits (argv[2]));
	    if (ERRPRNT ((fn_hr = (uint64_t) strtol (argv[2] + 10, NULL, 10),
			  errno != 0)))
	      {
		exit (1);
	      }
	  }
	else
	  {
	    fprintf (stderr, "\n*** invalid second argument ***\n\n");
	    fprintf (stderr, USAGE);
	    exit (1);
	  }
      }
      break;
    default:
      fprintf (stderr, "*** wrong number of arguments ***\n\n");
      fprintf (stderr, USAGE);
      exit (1);
    }
  /* END{process the command line and arguments} */

  /* now do the work */
  if (fn_hr > UCHAR_MAX)
    {
      FILE *fp;
      struct stat statbuf;
      set_hrbytes (fn_hr);
      ERREXIT (stat (fname, &statbuf) != 0);
      ERREXIT ((fp = fopen (fname, "r+b")) == NULL);
      while (! feof (fp)) {
  	if (gribmsgfixtime (fp) >= statbuf.st_size) break;
      }
      ERRPRNT (fclose (fp) != 0);
    }

  return 0;
}

static int
all_digits (char *s)
{
  int i;
  for (i = 0; s[i] != '\0'; i++)
    {
      if (! isdigit(s[i]))
	{
	  errno = EINVAL;
	  return 0;
	}
    }
  errno = 0;
  return 1;
}
