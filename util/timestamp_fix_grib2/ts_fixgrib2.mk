# -*- mode: makefile-gmake -*-

SHELL = /bin/sh

ts_fixgrib2_exe = timestamp_fix_grib2
ts_fixgrib2_src = timestamp_fix_grib2.c gribmsgfixtime.c
ts_fixgrib2_hdr = errmacs.h gribmsgfixtime.h
ts_fixgrib2_mk = ts_fixgrib2.mk
ts_fixgrib2_version = 2.0
ts_fixgrib2_obj = $(ts_fixgrib2_src:.c=.o)
ts_fixgrib2_dep = $(addprefix .,$(ts_fixgrib2_src:.c=.d))

# CDEPFLAG should be -MM for gcc, -M otherwise but still works with gcc
CDEPFLAG = -MM
# These following make variables are set with values that work with gcc
CFLAGS = -O3
CPPFLAGS = -Wall

all: $(ts_fixgrib2_exe)
$(ts_fixgrib2_exe): $(ts_fixgrib2_obj)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)
.%.d: %.c
	@set -e; rm -f $@; \
	$(CC) $(CDEPFLAG) $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
	@echo $*

tar.gz: timestamp_fix_grib2-$(ts_fixgrib2_version).tar.gz
timestamp_fix_grib2-$(ts_fixgrib2_version).tar.gz: $(ts_fixgrib2_src) \
                                                   $(ts_fixgrib2_hdr) \
                                                   $(ts_fixgrib2_mk)
	mkdir timestamp_fix_grib2-$(ts_fixgrib2_version)
	cp -p $^ timestamp_fix_grib2-$(ts_fixgrib2_version)/.
	cp -p $(ts_fixgrib2_mk) \
	      timestamp_fix_grib2-$(ts_fixgrib2_version)/Makefile
	tar covpf timestamp_fix_grib2-$(ts_fixgrib2_version).tar \
	          timestamp_fix_grib2-$(ts_fixgrib2_version)
	gzip -9 timestamp_fix_grib2-$(ts_fixgrib2_version).tar

.PHONY: clean
clean:
	$(RM) $(ts_fixgrib2_exe) $(ts_fixgrib2_obj) $(ts_fixgrib2_dep) *~
	$(RM) -r timestamp_fix_grib2-$(ts_fixgrib2_version)
	$(RM) timestamp_fix_grib2-$(ts_fixgrib2_version).tar.gz

include $(ts_fixgrib2_dep)
