#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "errmacs.h"

struct section_header_t {
  long int length;
  int number;
};
typedef struct section_header_t section_header_t;
static uint64_t forecast_time;
static uint8_t hrbyte1, hrbyte2;
static int gribmsg;
static long int toppos;
static long int msglen = 0L;

/* This takes up the first 5 bytes of a section. */
static section_header_t * section_header (FILE *fp, section_header_t *h);

void
set_hrbytes (uint64_t fn_hr)
{
  forecast_time = fn_hr;
  hrbyte1 = (uint8_t) (fn_hr >> 8);
  hrbyte2 = (uint8_t) (fn_hr & 0xFFL);
  gribmsg = 0;
}

#define go_past_sec(n)							\
  if (sec_hdr.number != n)						\
    {									\
      fprintf (stderr, "*** should be at section " #n " of record %d "	\
	       "not section %d.",					\
	       gribmsg, sec_hdr.number);				\
      fprintf (stderr, "  exiting...\n");				\
      exit (1);								\
    }									\
  /* already past first 5 (header) bytes */				\
  /* set file position to next section */				\
  err_fseek (fp, sec_hdr.length - 5, SEEK_CUR);				\
  /* read header for next section */					\
  section_header (fp, &sec_hdr);

long int
gribmsgfixtime (FILE *fp)
{
  section_header_t sec_hdr;
  uint8_t grib_ed;
  int i;
  uint8_t c;
  const uint8_t grib[] = "GRIB";
  uint8_t buf4ch[4];

  errno = 0;

  err_ftell (toppos, fp);
  gribmsg++;

  /* Section 0 */
  /* check for "GRIB" string at beginning of message */
  err_fread (buf4ch, sizeof (uint8_t), 4, fp);
  if (memcmp (grib, buf4ch, 4) != 0)
    {
      fprintf (stderr, "*** Can't find beginning of GRIB message\n");
      exit (1);
    }

  err_fseek (fp, 3, SEEK_CUR);	/* pass the eightth byte */
  err_getc (grib_ed, fp);
  /* fprintf (stderr, "GRIB edition %d\n", grib_ed); */
  /* get entire message length in number of bytes, bytes 9 to 16 */
  for (i = 0; i < 8; i++)
    {
      err_getc (c, fp);
      is_early_eof (c);
      msglen <<= CHAR_BIT;
      msglen |= (c & UCHAR_MAX);
    }

  section_header (fp, &sec_hdr);

  /* Section 1 */
  go_past_sec(1);

  while (1)
    {
      switch (sec_hdr.number)
	{
	case 2: /* Section 2 */
	  go_past_sec (2);

	case 3: /* Section 3 */
	  go_past_sec (3);

	case 4: /* Section 4 */
	  if (sec_hdr.number == 4)
	    {
	      uint16_t prod_def_tmpl;
	      /* already past first 5 (header) bytes */
	      /* number of coordinate values after template */
	      err_fseek (fp, 2, SEEK_CUR);
	      /* Product Definition Template number */
	      err_getc (c, fp);	/* octet 8 */
	      prod_def_tmpl = (c & UCHAR_MAX) << CHAR_BIT;
	      err_getc (c, fp);	/* octet 9 */
	      prod_def_tmpl |= (c & UCHAR_MAX);

	      if (prod_def_tmpl == 0) /* Prod Def Tmpl is 4.0 */
		{
		  err_fseek (fp, 11, SEEK_CUR);
		  err_getc (c, fp); /* octet 21 */
		  if ((c & UCHAR_MAX) == '\0')
		    {
		      err_getc (c, fp); /* octet 22 */
		      if ((c & UCHAR_MAX) == hrbyte2)
			{
			  err_fseek (fp, -2, SEEK_CUR);
			  err_putc (hrbyte1, fp); /* octet 21 */
			  /* set file position to next section */
			  err_fseek (fp, sec_hdr.length - 21, SEEK_CUR);
			}
		      else
			{
			  fprintf (stderr,
				   "WARNING --- octet 22 of rec %d was "
				   "expected to be %d but was %d\n",
				   gribmsg, hrbyte2, c);
			  /* set file position to next section */
			  err_fseek (fp, sec_hdr.length - 22, SEEK_CUR);
			}
		    }
		  else
		    {
		      err_fseek (fp, sec_hdr.length - 21, SEEK_CUR);
		    }
		} /* if (Prod Def Tmpl is 4.0) */

	      else if (prod_def_tmpl == 8) /* Prod Def Tmpl is 4.8 */
		{
		  uint64_t time_range;
		  uint64_t forecast_end;
		  uint64_t fixed_forecast_start;
		  uint64_t forecast_start = 0L;

		  err_fseek (fp, 9, SEEK_CUR); /* octet 19 */
		  for (i = 0; i < 4; i++)      /* read octets 19 through 22 */
		    {
		      forecast_start <<= CHAR_BIT;
		      err_getc (c, fp);
		      forecast_start |= (c & UCHAR_MAX);
		    }

		  /* Adjust forecast start by UCHAR_MAX, if necessary...
		     i.e., not zero, and not already adjusted */
		  if (forecast_start != 0)
		    { 
		      if (forecast_start < UCHAR_MAX - 24) forecast_start += UCHAR_MAX + 1;
		    }
		  /*fprintf (stderr, "final forecast_start is: %d\n", forecast_start); */

		  err_fseek (fp, 27, SEEK_CUR); /* octet 50 */
		  /* This next for loop overwrites the first 3 out of
		     4 bytes of the time range because they seem like
		     they always need to be zeroed out and the second
		     and third bytes seem to have their bits all set
		     to 1.  This step is questionable and will
		     probably have to addressed with a more acceptable
		     answer later. */
		  for (i = 0; i < 3; i++) /* write octets 50 through 53 */
		    {
		      err_putc ('\0', fp);
		    }
		  /* setting the time_range variable to only the
		     last of the 4 time range bytes, again
		     questionable */
		  err_getc (c, fp); /* read octet 53 */
		  time_range = 0xFFL & c;
		  /* Assuming forecast_start is zero or small enough */
		  if ((hrbyte1 > 0) &&
		      ((time_range | (hrbyte1 << CHAR_BIT)) == forecast_time))
		    {
		      time_range = forecast_time;
		      err_fseek (fp, -2, SEEK_CUR); /* go to octet 52 */
		      err_putc (hrbyte1, fp);	    /* write to octet 52 */
		      err_fseek (fp, 1, SEEK_CUR); /* at octet 53 after write;
						      go to octet 54 */
		    }
		  forecast_end = forecast_start + time_range;
		  fixed_forecast_start = forecast_time - time_range;
		  if (forecast_end == forecast_time)
		    {
		      /* octets 19 to 22 */
		      err_fseek (fp, -35, SEEK_CUR);
		      for (i = 3; i >= 0; i--)
			{
			  err_putc ((fixed_forecast_start >>
				     (CHAR_BIT * i)) & UCHAR_MAX, fp);
			}
		      /* set file position to next section */
		      err_fseek (fp, sec_hdr.length - 22, SEEK_CUR);
		    }
		  else
		    {
		      fprintf (stderr,
			       "WARNING --- discrepancy in rec %d:\n"
			       "forecast_end summed from forecast_start and "
			       "time_range is %ld\n"
			       "forecast_time obtained from file name is %ld\n"
			       "forecast_time and forecast_end should be "
			       "the same\n",
			       gribmsg, forecast_end, forecast_time);
		      /* set file position to next section */
		      err_fseek (fp, sec_hdr.length - 53, SEEK_CUR);
		    }
		} /* else if (Prod Def Tmpl is 4.8) */
	      else
		{
		  err_fseek (fp, sec_hdr.length - 9, SEEK_CUR);
		}
	      /* read header for the next section */
	      section_header (fp, &sec_hdr);
	    }

	  /* Section 5 */
	  go_past_sec (5);

	  /* Section 6 */
	  go_past_sec (6);

	  /* Section 7 */
	  go_past_sec (7);

	  /* Section 8 */
	  if (sec_hdr.number == 8)
	    {
	      err_ftell (toppos, fp);
	      return toppos;
	    }
	  break;		/* switch/case */

	default:
	  fprintf (stderr, "*** section number is %d\n", sec_hdr.number);
	  fprintf (stderr, "can't find proper section header.  exiting...\n");
	  exit (1);

	} /* switch (sec_hdr.number) */

    } /* while (1) */
}

long int
checkgribmsgtime (FILE *fp)
{
  long int curpos;
  long int msglen = 0L;
  int i;
  int c;
  const char grib[] = "GRIB";

  errno = 0;

  err_ftell (curpos, fp);
  gribmsg++;

  /* check for "GRIB" string at beginning of message */
  for (i = 0; i < 4; i++)
    {
      err_getc (c, fp);
      is_early_eof (c);
      if (grib[i] != c)
	{
	  fprintf (stderr, "*** Can't find beginning of GRIB message\n");
	  exit (1);
	}
    }

  /* get entire message length in number of bytes, bytes 9 to 16 */
  err_fseek (fp, 4, SEEK_CUR);	/* pass the eightth byte */
  for (i = 0; i < 8; i++)
    {
      err_getc (c, fp);
      is_early_eof (c);
      msglen <<= CHAR_BIT;
      msglen |= (c & UCHAR_MAX);
    }

  /* check for "7777" string at end of message */
  err_fseek (fp, msglen - 20, SEEK_CUR); /* go back for the trailing
					    "7777" */
  for (i = 0; i < 4; i++)
    {
      err_getc (c, fp);
      if ('7' != c)
	{
	  fprintf (stderr, "%d\n", i);
	  fprintf (stderr, "*** Can't find end of GRIB message\n");
	  exit (1);
	}
    }

  err_ftell (curpos, fp);
  return curpos;
}

/* This takes up the first 5 bytes of a section. */
static section_header_t *
section_header (FILE *fp, section_header_t *h)
{
  char clenbuf[4];
  long int curpos;
  int i;
  err_fread (clenbuf, sizeof (char), 4, fp);
  if (strncmp (clenbuf, "7777", 4) == 0)
    {
      h->number = 8;
      h->length = 4L;
      return h;
    }
  h->length = 0L;
  for (i = 0; i < 4; i++)
    {
      h->length <<= CHAR_BIT;
      h->length |= (clenbuf[i] & UCHAR_MAX);
    }
  err_getc (h->number, fp);
  err_ftell (curpos, fp);
  if (curpos - toppos > msglen)
    {
      fprintf (stderr, "read past end of message.  exiting...\n");
      exit (1);
    }
  return h;
}
