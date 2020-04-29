#ifndef ERRMACS_H
#define ERRMACS_H 1

#include <stdlib.h>
#include <errno.h>

#define ERREXIT(X)							\
  if ((X) != 0) {							\
    fprintf (stderr, "Near line %d of file %s\n", __LINE__, __FILE__);	\
    perror(#X);								\
    exit (1);								\
  }
#define ERRPRNT(X) ((X) != 0 ? (perror(#X), errno = 0, 1) : 0)

#define gribmsg_ERRVAL -1L
#define handle_ferror(fp)						\
  if (ferror (fp))							\
    {									\
      fprintf (stderr, "*** Error on file stream.  Exiting...\n");	\
      exit (1);								\
    }
#define is_early_eof(c)							\
  if ((c) == EOF)							\
    {									\
      fprintf (stderr, "WARNING - end of file encountered...\n");	\
      fprintf (stderr, "   ---   near line %d of file %s.\n",		\
	       __LINE__, __FILE__);					\
      return gribmsg_ERRVAL;						\
    }
#define err_getc(c,fp)				\
  (c) = getc (fp);				\
  if ((c) == EOF)				\
    {						\
      handle_ferror (fp);			\
    }
#define err_fread(ptr,sz,num,fp)					\
  fread ((ptr), (sz), (num), (fp));					\
  if (ferror (fp))							\
    {									\
      fprintf (stderr, "*** Error on file stream.  Exiting...\n");	\
      exit (1);								\
    }									\
  else if (feof (fp))							\
    {									\
      fprintf (stderr, "WARNING - end of file encountered...\n");	\
      fprintf (stderr, "   ---   near line %d of file %s.\n",		\
	       __LINE__, __FILE__);					\
    }
#define err_ungetc(c,fp)			\
  ungetc ((c), (fp));				\
  handle_ferror (fp)
#define err_putc(c,fp)				\
  if (putc ((c), (fp)) == EOF)			\
    {						\
      handle_ferror (fp);			\
    }
#define err_fseek(fp,o,w)			\
  fseek ((fp), (o), (w));			\
  handle_ferror (fp)
#define err_ftell(p,fp)				\
  (p) = ftell (fp);				\
  handle_ferror (fp)

#endif	/* ERRMACS_H */
