/*
** svn $Id: pete_NWP6.h 139 2008-01-10 00:17:29Z arango $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Norh West Pacific case.
**
** Application flag:   ecsy10 from NWP4_level 10
** Input script:       ecsy10.in
*/

#define WET_DRY
#define UV_ADV
#define UV_QDRAG
#define UV_COR
#define UV_VIS2 
#define TS_DIF2 
#define MIX_S_UV
#define MIX_GEO_TS
#define DJ_GRADPS
#define TS_U3HADVECTION
#define NONLIN_EOS
#define SALINITY
#define SOLVE3D
#define MASKING
#define SPLINES_VVISC
#define SPLINES_VDIFF
#define ANA_SPONGE

#undef QCORRECTION  /*net heat flux correction */
#define SCORRECTION  /*freshwater flux correction*/
#define SOLAR_SOURCE /*solar radiation source term*/
#undef DIURNAL_SRFLUX /*mdulate input shortwave with diurnal cycle*/
#undef SRELAXATION   /*salinity relaxation as a freshwater flux*/
/*BULK FLUX*/

#define BULK_FLUXES
#define LONGWAVE_OUT
#define EMINUSP
#define LIMIT_STFLX_COOLING

#define CURVGRID
#undef DIAGNOSTICS_TS
#undef DIAGNOSTICS_UV
#define INLINE_2DIO
#undef ATM_PRESS
/* output*/

#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS
#define AVERAGES_FLUXES

/*#define AVERAGES_DETIDE*/
#define AVERAGES
#undef STATIONS
#undef FLOATS
#undef PERFECT_RESTART

/* option for vertical mixing */
#undef GLS_MIXING
#undef MY25_MIXING
#define LMD_MIXING

#ifdef GLS_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
# define RI_SPLINES
#endif

#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
#endif

#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# undef LMD_DDMIX
# define LMD_SKPP
# undef LMD_BKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif

/* set bottom boundary value zero*/
#define ANA_BSFLUX
#define ANA_BTFLUX


#undef UV_TIDES
#undef SSH_TIDES
#undef ADD_M2OBC
#undef ADD_FSOBC
#undef RAMP_TIDES

#define POSITIVEDEF
#undef NPZD_IRON_RESTRUCT
#define NPZD_IRON
# ifdef NPZD_IRON
#  define IRON_LIMIT
#  define IRON_RELAX
# endif
#define ANA_SPFLUX
#define ANA_BPFLUX

