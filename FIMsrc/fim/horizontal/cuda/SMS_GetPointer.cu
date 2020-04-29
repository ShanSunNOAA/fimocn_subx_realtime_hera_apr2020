// This routine is used by SMS.
// This routine returns a pointer to an exchange variable.
// The pointers table in this routine must contain all the variables to be exchanged.
// Eventually PPP will generate the pointers table.
// Currently (April 2012) the pointers table is hard coded for FIM.
// Author:  Jacques Middlecoff
// Date:    August 2012

// Variables that appear here must be malloc'd in one of the copytoGPU routines.

#include <stdio.h>
#include <string.h>
extern "C" void sms_getpointer (char *varname, void **ptr, int *status, int string_len) {

extern int   *D_prox;
extern float *D_nprox;
extern float *D_proxs;
extern float *D_lat;
extern float *D_lon;
extern float *D_delp_lo;
extern float *D_delp;
extern float *D_r_plus;
extern float *D_r_mnus;
extern float *D_ptdcy;
extern float *D_area;
extern float *D_cs;
extern float *D_sn;
extern float *D_sidevec_c;
extern float *D_sidevec_e;
extern float *D_sideln;
extern float *D_rprox_ln;
extern float *D_corio;
extern float *D_deg_lat;
extern float *D_deg_lon;
extern float *D_rarea;
extern float *D_rsideln;
extern float *D_actual;
extern float *D_work_edg;
extern float *D_us3d;
extern float *D_vs3d;
extern float *D_dp3d;
extern float *D_tr3d;
extern float *D_mp3d;
extern float *D_pr3d;
extern float *D_fields;
extern float *D_thko;
extern float *D_thkn;
extern float *D_flp;
extern float *D_fln;
extern float *D_workb;
extern float *D_exsmo3d;
extern float *D_fld;
extern float *D_cumufx;
extern float *D_dpinit;
extern int   *D_dpfinl;
extern float *D_massfx;
extern float *D_trcr_edg;
extern float *D_tracr;
extern float *D_trcr_lo;
extern float *D_g3p;

struct pointers {
  char *string;
  void *D_;
};

struct pointers table[] = {{"prox"     , D_prox     },
                           {"nprox"    , D_nprox    },
                           {"proxs"    , D_proxs    },
                           {"lat"      , D_lat      },
                           {"lon"      , D_lon      },
                           {"delp_lo"  , D_delp_lo  },
                           {"delp"     , D_delp     },
                           {"r_plus"   , D_r_plus   },
                           {"r_mnus"   , D_r_mnus   },
                           {"ptdcy"    , D_ptdcy    },
                           {"area"     , D_area     },
                           {"cs"       , D_cs       },
                           {"sn"       , D_sn       },
                           {"sidevec_c", D_sidevec_c},
                           {"sidevec_e", D_sidevec_e},
                           {"sideln"   , D_sideln   },
                           {"rprox_ln" , D_rprox_ln },
                           {"corio"    , D_corio    },
                           {"deg_lat"  , D_deg_lat  },
                           {"deg_lon"  , D_deg_lon  },
                           {"rarea"    , D_rarea    },
                           {"rsideln"  , D_rsideln  },
                           {"actual"   , D_actual   },
                           {"work_edg" , D_work_edg },
                           {"us3d"     , D_us3d     },
                           {"vs3d"     , D_vs3d     },
                           {"dp3d"     , D_dp3d     },
                           {"tr3d"     , D_tr3d     },
                           {"mp3d"     , D_mp3d     },
                           {"pr3d"     , D_pr3d     },
                           {"fields"   , D_fields   },
                           {"thko"     , D_thko     },
                           {"thkn"     , D_thkn     },
                           {"flp"      , D_flp      },
                           {"fln"      , D_fln      },
                           {"workb"    , D_workb    },
                           {"exsmo3d"  , D_exsmo3d  },
                           {"fld"      , D_fld      },
                           {"cumufx"   , D_cumufx   },
                           {"dpinit"   , D_dpinit   },
                           {"dpfinl"   , D_dpfinl   },
                           {"massfx"   , D_massfx   },
                           {"trcr_edg" , D_trcr_edg },
                           {"tracr"    , D_tracr    },
                           {"trcr_lo"  , D_trcr_lo  },
                           {"g3p"      , D_g3p      }
                          };

 int i;
 const int numelem = sizeof (table) / sizeof (struct pointers);

 *status=0;

  for (i = 0; i < numelem; ++i) {
//  printf("JFM: %d %s %d %d %s %d \n",i,varname,string_len,strlen(table[i].string),table[i].string,status);
    if( (string_len == strlen(table[i].string)) && (strncmp(varname,table[i].string,string_len) == 0) ) { 
      *ptr = table[i].D_ ;
      return;
    }
  }
  printf("Error in SMS_GetPointer.cu: VarName out of range %s %d \n",varname,string_len);
  printf("Make sure all exchange variables are in the pointers table in SMS_GetPointer.cu \n");
  *status=-8003;
  return;
}
