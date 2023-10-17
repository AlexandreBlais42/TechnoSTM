/*
 *  $Id: gwyshapefitpreset.c 24952 2022-08-26 15:24:56Z yeti-dn $
 *  Copyright (C) 2016-2022 David Necas (Yeti).
 *  E-mail: yeti@gwyddion.net.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 *  License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any
 *  later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not, write to the
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "config.h"
#include <string.h>
#include <libgwyddion/gwymacros.h>
#include <libgwyddion/gwymath.h>
#include <libgwyddion/gwyrandgenset.h>
#include <libgwyddion/gwydebugobjects.h>
#include <libprocess/dataline.h>
#include <libprocess/linestats.h>
#include <libprocess/surface.h>
#include <libprocess/peaks.h>
#include <libprocess/gwyshapefitpreset.h>

enum { NREDLIM = 4096 };

/* Lower symmetric part indexing */
/* i MUST be greater or equal than j */
#define SLi(a, i, j) a[(i)*((i) + 1)/2 + (j)]

#define DECLARE_SECONDARY(funcname,name) \
    static gdouble funcname##_calc_##name    (const gdouble *param); \
    static gdouble funcname##_calc_err_##name(const gdouble *param, \
                                              const gdouble *param_err, \
                                              const gdouble *correl);

#define DECLARE_SHAPE_FUNC(name) \
    static gdouble name##_func(gdouble x, \
                               gdouble y, \
                               const gdouble *param); \
    static gdouble name##_fitfunc(guint i, \
                                  const gdouble *param, \
                                  gpointer user_data, \
                                  gboolean *fres) \
    { \
        const GwyXYZ *xyz = ((const ShapeFitPreset*)user_data)->xyz; \
        *fres = TRUE; \
        return name##_func(xyz[i].x, xyz[i].y, param) - xyz[i].z; \
    } \
    static gboolean name##_estimate(const GwyXYZ *xyz, \
                                    guint n, \
                                    gdouble *param, \
                                    FitShapeEstimateCache *estimcache); \
    static gboolean name##_init(const GwyXYZ *xyz, \
                                guint n, \
                                gdouble *param, \
                                FitShapeEstimateCache *estimcache);

#define DECLARE_SHAPE_DIFF(name) \
    static void name##_diff(gdouble x, \
                            gdouble y, \
                            const gdouble *param, \
                            const gboolean *fixed_param, \
                            gdouble *der); \
    static void name##_difffunc(guint i, \
                                const gdouble *param, \
                                const gboolean *fixed_param, \
                                G_GNUC_UNUSED GwyNLFitIdxFunc func, \
                                gpointer user_data, \
                                gdouble *der, \
                                gboolean *fres) \
    { \
        const GwyXYZ *xyz = ((const ShapeFitPreset*)user_data)->xyz; \
        *fres = TRUE; \
        name##_diff(xyz[i].x, xyz[i].y, param, fixed_param, der); \
    }

/* XXX: This is a dirty trick assuming sizeof(FitShapeSecondary) > sizeof(NULL) so that we get zero nsecondary when
 * name##_secondary is defined to NULL and correct array size otherwise.  It is safe because FitShapeSecondary is
 * a struct that contains at least two pointers plus other stuff, but it is dirty anyway.  Static code analysis tends
 * to dislike it. */
#define SHAPE_FUNC_ITEM(name) \
    &name##_func, &name##_fitfunc, NULL, &name##_estimate, &name##_init, \
    G_N_ELEMENTS(name##_params), sizeof(name##_secondary)/sizeof(FitShapeSecondary), name##_params, \
    name##_secondary

#define SHAPE_FDIF_ITEM(name) \
    &name##_func, &name##_fitfunc, &name##_difffunc, &name##_estimate, &name##_init, \
    G_N_ELEMENTS(name##_params), sizeof(name##_secondary)/sizeof(FitShapeSecondary), name##_params, \
    name##_secondary

typedef struct {
    gboolean have_mean;
    gboolean have_circle;
    gboolean have_zrange;
    gboolean have_zstats;
    /* Plain mean values */
    gdouble xm;
    gdouble ym;
    /* Circumscribed circle. */
    gdouble xc;
    gdouble yc;
    gdouble r;
    /* Value range. */
    gdouble zmin;
    gdouble zmax;
    /* Simple value stats. */
    gdouble zmean;
    gdouble zrms;
    gdouble zskew;
} FitShapeEstimateCache;

typedef gdouble (*FitShapeXYFunc)(gdouble x, gdouble y,
                                  const gdouble *param);

typedef gboolean (*FitShapeEstimate)(const GwyXYZ *xyz,
                                     guint n,
                                     gdouble *param,
                                     FitShapeEstimateCache *estimcache);

typedef gdouble (*FitShapeCalcParam)(const gdouble *param);
typedef gdouble (*FitShapeCalcError)(const gdouble *param,
                                     const gdouble *param_err,
                                     const gdouble *correl);

typedef struct {
    const char *name;
    const char *desc;
    gint power_xy;
    gint power_z;
    GwyNLFitParamFlags flags;
} FitShapeParam;

typedef struct {
    const char *name;
    const char *desc;
    gint power_xy;
    gint power_z;
    GwyNLFitParamFlags flags;
    FitShapeCalcParam calc;
    FitShapeCalcError calc_err;
} FitShapeSecondary;

typedef struct {
    const gchar *name;
    gboolean needs_same_units;
    FitShapeXYFunc function;
    GwyNLFitIdxFunc fit_function;
    GwyNLFitIdxDiffFunc diff_function;
    FitShapeEstimate estimate;
    FitShapeEstimate initialise;
    guint nparams;
    guint nsecondary;
    const FitShapeParam *param;
    const FitShapeSecondary *secondary;
} FitShapeFunc;

struct _GwyShapeFitPresetPrivate {
    const FitShapeFunc *builtin;
    /* Context for the opaque function during fitting. */
    const GwyXYZ *xyz;
};

typedef struct _GwyShapeFitPresetPrivate ShapeFitPreset;

DECLARE_SHAPE_FUNC(plane);
DECLARE_SHAPE_DIFF(plane);
DECLARE_SHAPE_FUNC(step1);
DECLARE_SHAPE_DIFF(step1);
DECLARE_SHAPE_FUNC(step2);
DECLARE_SHAPE_FUNC(grating);
DECLARE_SHAPE_FUNC(grating3);
DECLARE_SHAPE_FUNC(holes);
DECLARE_SHAPE_FUNC(sphere);
DECLARE_SHAPE_FUNC(hsphere);
DECLARE_SHAPE_FUNC(hcylinder);
DECLARE_SHAPE_FUNC(gaussian);
DECLARE_SHAPE_DIFF(gaussian);
DECLARE_SHAPE_FUNC(lorentzian);
DECLARE_SHAPE_DIFF(lorentzian);
DECLARE_SHAPE_FUNC(exponential);
DECLARE_SHAPE_DIFF(exponential);
DECLARE_SHAPE_FUNC(pyramid);
DECLARE_SHAPE_FUNC(pyramidx);
DECLARE_SHAPE_FUNC(pyramid3);
DECLARE_SHAPE_FUNC(parbump);
DECLARE_SHAPE_FUNC(cone);
DECLARE_SHAPE_FUNC(smthcone);
DECLARE_SHAPE_FUNC(smthpyr3);
DECLARE_SHAPE_FUNC(smthpyr4);
DECLARE_SHAPE_FUNC(cylinder);
DECLARE_SHAPE_DIFF(cylinder);
DECLARE_SHAPE_FUNC(pring);
DECLARE_SHAPE_FUNC(lunette);

DECLARE_SECONDARY(plane, phi);
DECLARE_SECONDARY(plane, theta);
DECLARE_SECONDARY(grating3, h);
DECLARE_SECONDARY(grating3, L0);
DECLARE_SECONDARY(grating3, L1);
DECLARE_SECONDARY(grating3, L2);
DECLARE_SECONDARY(grating3, L3);
DECLARE_SECONDARY(holes, wouter);
DECLARE_SECONDARY(holes, winner);
DECLARE_SECONDARY(holes, R);
DECLARE_SECONDARY(sphere, R);
DECLARE_SECONDARY(sphere, zc);
DECLARE_SECONDARY(hsphere, R);
DECLARE_SECONDARY(hsphere, zc);
DECLARE_SECONDARY(hcylinder, R);
DECLARE_SECONDARY(gaussian, sigma1);
DECLARE_SECONDARY(gaussian, sigma2);
DECLARE_SECONDARY(lorentzian, b1);
DECLARE_SECONDARY(lorentzian, b2);
DECLARE_SECONDARY(exponential, b1);
DECLARE_SECONDARY(exponential, b2);
DECLARE_SECONDARY(cylinder, R1);
DECLARE_SECONDARY(cylinder, R2);

static const FitShapeParam plane_params[] = {
   { "z<sub>0</sub>", N_("Value at the origin"),      0,  1, 0, },
   { "b<sub>x</sub>", N_("Base plane x-coefficient"), -1, 1, 0, },
   { "b<sub>y</sub>", N_("Base plane y-coefficient"), -1, 1, 0, },
};

static const FitShapeSecondary plane_secondary[] = {
   { "φ", N_("Angle from horizontal plane"), 0, 0, GWY_NLFIT_PARAM_ANGLE, plane_calc_phi,   plane_calc_err_phi,   },
   { "ϑ", N_("Orientation"),                 0, 0, GWY_NLFIT_PARAM_ANGLE, plane_calc_theta, plane_calc_err_theta, },
};

static const FitShapeParam step1_params[] = {
   { "h",             N_("Step height"),                             0,  1, 0,                      },
   { "δ",             N_("Slope width"),                             1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "z<sub>0</sub>", N_("Base plane value at the origin"),          0,  1, 0,                      },
   { "x<sub>0</sub>", N_("Distance of step center from the origin"), 1,  0, 0,                      },
   { "φ",             N_("Direction perpendicular to edges"),        0,  0, GWY_NLFIT_PARAM_ANGLE,  },
   { "b<sub>x</sub>", N_("Base plane x-coefficient"),                -1, 1, 0,                      },
   { "b<sub>y</sub>", N_("Base plane y-coefficient"),                -1, 1, 0,                      },
};

#define step1_secondary NULL

static const FitShapeParam step2_params[] = {
   { "h",             N_("Step height"),                             0,  1, 0,                      },
   { "w",             N_("Width at half-height"),                    1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "δ",             N_("Slope width"),                             1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "z<sub>0</sub>", N_("Base plane value at the origin"),          0,  1, 0,                      },
   { "x<sub>0</sub>", N_("Distance of step center from the origin"), 1,  0, 0,                      },
   { "φ",             N_("Direction perpendicular to edges"),        0,  0, GWY_NLFIT_PARAM_ANGLE,  },
   { "b<sub>x</sub>", N_("Base plane x-coefficient"),                -1, 1, 0,                      },
   { "b<sub>y</sub>", N_("Base plane y-coefficient"),                -1, 1, 0,                      },
};

#define step2_secondary NULL

static const FitShapeParam grating_params[] = {
   { "L",             N_("Spatial period (pitch)"),           1, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "h",             N_("Step height"),                      0, 1, 0,                      },
   { "p",             N_("Fill ratio at base"),               0, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "z<sub>0</sub>", N_("Base plane value at the origin"),   0, 1, 0,                      },
   { "x<sub>0</sub>", N_("Offset from origin"),               1, 0, 0,                      },
   { "φ",             N_("Direction perpendicular to edges"), 0, 0, GWY_NLFIT_PARAM_ANGLE,  },
   { "c",             N_("Top edge sharpness"),               0, 0, GWY_NLFIT_PARAM_ABSVAL, },
};

#define grating_secondary NULL

static const FitShapeParam grating3_params[] = {
   { "L",             N_("Spatial period (pitch)"),           1, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "h<sub>1</sub>", N_("Bottom tier height"),               0, 1, GWY_NLFIT_PARAM_ABSVAL, },
   { "h<sub>2</sub>", N_("Middle tier height"),               0, 1, GWY_NLFIT_PARAM_ABSVAL, },
   { "h<sub>3</sub>", N_("Top tier height"),                  0, 1, GWY_NLFIT_PARAM_ABSVAL, },
   { "p",             N_("Fill ratio at base"),               0, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "q<sub>1</sub>", N_("Bottom tier tapering factor"),      0, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "q<sub>2</sub>", N_("Middle tier tapering factor"),      0, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "q<sub>3</sub>", N_("Top tier tapering factor"),         0, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "z<sub>0</sub>", N_("Base plane value at the origin"),   0, 1, 0,                      },
   { "x<sub>0</sub>", N_("Offset from origin"),               1, 0, 0,                      },
   { "φ",             N_("Direction perpendicular to edges"), 0, 0, GWY_NLFIT_PARAM_ANGLE,  },
};

static const FitShapeSecondary grating3_secondary[] = {
   { "h",             N_("Feature height"), 0, 1, 0, grating3_calc_h,  grating3_calc_err_h,  },
   { "L<sub>0</sub>", NULL,                 1, 0, 0, grating3_calc_L0, grating3_calc_err_L0, },
   { "L<sub>1</sub>", NULL,                 1, 0, 0, grating3_calc_L1, grating3_calc_err_L1, },
   { "L<sub>2</sub>", NULL,                 1, 0, 0, grating3_calc_L2, grating3_calc_err_L2, },
   { "L<sub>3</sub>", NULL,                 1, 0, 0, grating3_calc_L3, grating3_calc_err_L3, },
};

static const FitShapeParam holes_params[] = {
   { "x<sub>0</sub>", N_("Offset from origin in X"),        1, 0, 0,                      },
   { "y<sub>0</sub>", N_("Offset from origin in Y"),        1, 0, 0,                      },
   { "z<sub>0</sub>", N_("Base plane value at the origin"), 0, 1, 0,                      },
   { "L",             N_("Spatial period (pitch)"),         1, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "p",             NULL,                                 0, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "h",             N_("Step height"),                    0, 1, 0,                      },
   { "s",             NULL,                                 0, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "r",             N_("Corner roundness"),               0, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "φ",             N_("Orientation"),                    0, 0, GWY_NLFIT_PARAM_ANGLE,  },
};

static const FitShapeSecondary holes_secondary[] = {
   { "w<sub>outer</sub>", N_("Width at base"),    1, 0, 0, holes_calc_wouter, holes_calc_err_wouter, },
   { "w<sub>inner</sub>", N_("Width at the top"), 1, 0, 0, holes_calc_winner, holes_calc_err_winner, },
   { "R",                 N_("Corner radius"),    1, 0, 0, holes_calc_R,      holes_calc_err_R,      },
};

static const FitShapeParam sphere_params[] = {
   { "x<sub>0</sub>", N_("X-coordinate of the center"),     1, 0,  0, },
   { "y<sub>0</sub>", N_("Y-coordinate of the center"),     1, 0,  0, },
   { "z<sub>0</sub>", N_("Base plane value at the origin"), 0, 1,  0, },
   { "C",             N_("Curvature"),                      0, -1, 0, },
};

static const FitShapeSecondary sphere_secondary[] = {
   { "R",             N_("Radius"),             0, 1, 0, sphere_calc_R,  sphere_calc_err_R,  },
   { "z<sub>c</sub>", N_("Value at the apex"), 0, 1, 0, sphere_calc_zc, sphere_calc_err_zc, },
};

static const FitShapeParam hsphere_params[] = {
   { "x<sub>0</sub>", N_("X-coordinate of the center"),     1, 0,  0, },
   { "y<sub>0</sub>", N_("Y-coordinate of the center"),     1, 0,  0, },
   { "z<sub>0</sub>", N_("Base plane value at the origin"), 0, 1,  0, },
   { "C",             N_("Curvature"),                      0, -1, 0, },
};

static const FitShapeSecondary hsphere_secondary[] = {
   { "R",             NULL, 0, 1, 0, hsphere_calc_R,  hsphere_calc_err_R,  },
   { "z<sub>c</sub>", NULL, 0, 1, 0, hsphere_calc_zc, hsphere_calc_err_zc, },
};

static const FitShapeParam hcylinder_params[] = {
   { "x<sub>0</sub>", N_("Offset from origin"),                1,  0,  0,                     },
   { "z<sub>0</sub>", N_("Base plane value at the origin"),    0,  1,  0,                     },
   { "C",             N_("Curvature"),                         0,  -1, 0,                     },
   { "φ",             N_("Orientation"),                       0,  0,  GWY_NLFIT_PARAM_ANGLE, },
   { "b<sub>∥</sub>", N_("Base plane longitudal coefficient"), -1, 1,  0,                     },
};

static const FitShapeSecondary hcylinder_secondary[] = {
   { "R", NULL, 0, 1, 0, hcylinder_calc_R, hcylinder_calc_err_R, },
};

static const FitShapeParam gaussian_params[] = {
   { "x<sub>0</sub>",    N_("X-coordinate of the center"),                 1,  0, 0,                      },
   { "y<sub>0</sub>",    N_("Y-coordinate of the center"),                 1,  0, 0,                      },
   { "z<sub>0</sub>",    N_("Base plane value at the origin"),             0,  1, 0,                      },
   { "h",                N_("Height of the feature"),                      0,  1, 0,                      },
   { "σ<sub>mean</sub>", N_("Mean standard deviation"),                    1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "a",                N_("Elongation ratio"),                           0,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "φ",                N_("Direction of elongation (φ=0 means y-axis)"), 0,  0, GWY_NLFIT_PARAM_ANGLE,  },
   { "b<sub>x</sub>",    N_("Base plane x-coefficient"),                   -1, 1, 0,                      },
   { "b<sub>y</sub>",    N_("Base plane y-coefficient"),                   -1, 1, 0,                      },
};

static const FitShapeSecondary gaussian_secondary[] = {
   { "σ<sub>1</sub>", N_("Standard deviation 1"), 1, 0, 0, gaussian_calc_sigma1, gaussian_calc_err_sigma1, },
   { "σ<sub>2</sub>", N_("Standard deviation 2"), 1, 0, 0, gaussian_calc_sigma2, gaussian_calc_err_sigma2, },
};

static const FitShapeParam lorentzian_params[] = {
   { "x<sub>0</sub>",    N_("X-coordinate of the center"),                 1,  0, 0,                      },
   { "y<sub>0</sub>",    N_("Y-coordinate of the center"),                 1,  0, 0,                      },
   { "z<sub>0</sub>",    N_("Base plane value at the origin"),             0,  1, 0,                      },
   { "h",                N_("Height of the feature"),                      0,  1, 0,                      },
   { "β<sub>mean</sub>", N_("Mean width"),                                 1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "a",                N_("Elongation ratio"),                           0,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "φ",                N_("Direction of elongation (φ=0 means y-axis)"), 0,  0, GWY_NLFIT_PARAM_ANGLE,  },
   { "b<sub>x</sub>",    N_("Base plane x-coefficient"),                   -1, 1, 0,                      },
   { "b<sub>y</sub>",    N_("Base plane y-coefficient"),                   -1, 1, 0,                      },
};

static const FitShapeSecondary lorentzian_secondary[] = {
   { "β<sub>1</sub>", N_("Width 1"), 1, 0, 0, lorentzian_calc_b1, lorentzian_calc_err_b1, },
   { "β<sub>2</sub>", N_("Width 2"), 1, 0, 0, lorentzian_calc_b2, lorentzian_calc_err_b2, },
};

#define exponential_params lorentzian_params

static const FitShapeSecondary exponential_secondary[] = {
   { "β<sub>1</sub>", N_("Width 1"), 1, 0, 0, exponential_calc_b1, exponential_calc_err_b1, },
   { "β<sub>2</sub>", N_("Width 2"), 1, 0, 0, exponential_calc_b2, exponential_calc_err_b2, },
};

static const FitShapeParam pyramid_params[] = {
   { "x<sub>0</sub>", N_("X-coordinate of the center"),                 1,  0, 0,                      },
   { "y<sub>0</sub>", N_("Y-coordinate of the center"),                 1,  0, 0,                      },
   { "z<sub>0</sub>", N_("Base plane value at the origin"),             0,  1, 0,                      },
   { "h",             N_("Height of the feature"),                      0,  1, 0,                      },
   { "L",             N_("Side length of the base"),                    1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "a",             N_("Elongation ratio"),                           0,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "φ",             N_("Direction of elongation (φ=0 means y-axis)"), 0,  0, GWY_NLFIT_PARAM_ANGLE,  },
   { "b<sub>x</sub>", N_("Base plane x-coefficient"),                   -1, 1, 0,                      },
   { "b<sub>y</sub>", N_("Base plane y-coefficient"),                   -1, 1, 0,                      },
};

#define pyramid_secondary NULL

#define pyramidx_params pyramid_params
#define pyramidx_secondary NULL

#define pyramid3_params pyramid_params
#define pyramid3_secondary NULL

static const FitShapeParam cone_params[] = {
   { "x<sub>0</sub>", N_("X-coordinate of the center"),                 1,  0, 0,                      },
   { "y<sub>0</sub>", N_("Y-coordinate of the center"),                 1,  0, 0,                      },
   { "z<sub>0</sub>", N_("Base plane value at the origin"),             0,  1, 0,                      },
   { "h",             N_("Height of the feature"),                      0,  1, 0,                      },
   { "R",             N_("Base radius"),                                1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "a",             N_("Elongation ratio"),                           0,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "φ",             N_("Direction of elongation (φ=0 means y-axis)"), 0,  0, GWY_NLFIT_PARAM_ANGLE,  },
   { "b<sub>x</sub>", N_("Base plane x-coefficient"),                   -1, 1, 0,                      },
   { "b<sub>y</sub>", N_("Base plane y-coefficient"),                   -1, 1, 0,                      },
};

#define cone_secondary NULL

static const FitShapeParam smthcone_params[] = {
   { "x<sub>0</sub>", N_("X-coordinate of the center"),     1,  0, 0,                      },
   { "y<sub>0</sub>", N_("Y-coordinate of the center"),     1,  0, 0,                      },
   { "z<sub>0</sub>", N_("Base plane value at the origin"), 0,  1, 0,                      },
   { "h",             N_("Height of the feature"),          0,  1, 0,                      },
   { "α",             N_("Slope"),                          -1, 1, GWY_NLFIT_PARAM_ABSVAL, },
   { "R",             N_("Radius at apex"),                 1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "b<sub>x</sub>", N_("Base plane x-coefficient"),       -1, 1, 0,                      },
   { "b<sub>y</sub>", N_("Base plane y-coefficient"),       -1, 1, 0,                      },
};

#define smthcone_secondary NULL

static const FitShapeParam smthpyr3_params[] = {
   { "x<sub>0</sub>", N_("X-coordinate of the center"),                 1,  0, 0,                      },
   { "y<sub>0</sub>", N_("Y-coordinate of the center"),                 1,  0, 0,                      },
   { "z<sub>0</sub>", N_("Base plane value at the origin"),             0,  1, 0,                      },
   { "h",             N_("Height of the feature"),                      0,  1, 0,                      },
   { "α",             N_("Slope"),                                      -1, 1, GWY_NLFIT_PARAM_ABSVAL, },
   { "R",             N_("Radius at apex"),                             1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "φ",             N_("Direction of elongation (φ=0 means y-axis)"), 0,  0, GWY_NLFIT_PARAM_ANGLE,  },
   { "b<sub>x</sub>", N_("Base plane x-coefficient"),                   -1, 1, 0,                      },
   { "b<sub>y</sub>", N_("Base plane y-coefficient"),                   -1, 1, 0,                      },
};

#define smthpyr3_secondary NULL

#define smthpyr4_params smthpyr3_params
#define smthpyr4_secondary NULL

static const FitShapeParam parbump_params[] = {
   { "x<sub>0</sub>",    N_("X-coordinate of the center"),                 1,  0, 0,                      },
   { "y<sub>0</sub>",    N_("Y-coordinate of the center"),                 1,  0, 0,                      },
   { "z<sub>0</sub>",    N_("Base plane value at the origin"),             0,  1, 0,                      },
   { "h",                N_("Height of the feature"),                      0,  1, 0,                      },
   { "C<sub>mean</sub>", N_("Mean curvature at apex"),                     -1, 0, GWY_NLFIT_PARAM_ABSVAL, },
   { "a",                N_("Elongation ratio"),                           0,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "φ",                N_("Direction of elongation (φ=0 means y-axis)"), 0,  0, GWY_NLFIT_PARAM_ANGLE,  },
   { "b<sub>x</sub>",    N_("Base plane x-coefficient"),                   -1, 1, 0,                      },
   { "b<sub>y</sub>",    N_("Base plane y-coefficient"),                   -1, 1, 0,                      },
};

#define parbump_secondary NULL

static const FitShapeParam cylinder_params[] = {
   { "x<sub>0</sub>", N_("X-coordinate of the center"),                 1,  0, 0,                      },
   { "y<sub>0</sub>", N_("Y-coordinate of the center"),                 1,  0, 0,                      },
   { "z<sub>0</sub>", N_("Base plane value at the origin"),             0,  1, 0,                      },
   { "h",             N_("Height of the feature"),                      0,  1, 0,                      },
   { "R",             N_("Mean radius"),                                1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "C",             NULL,                                             0,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "a",             N_("Elongation ratio"),                           0,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "φ",             N_("Direction of elongation (φ=0 means y-axis)"), 0,  0, GWY_NLFIT_PARAM_ANGLE,  },
   { "b<sub>x</sub>", N_("Base plane x-coefficient"),                   -1, 1, 0,                      },
   { "b<sub>y</sub>", N_("Base plane y-coefficient"),                   -1, 1, 0,                      },
};

static const FitShapeSecondary cylinder_secondary[] = {
   { "R<sub>1</sub>", N_("Radius 1"), 1, 0, 0, cylinder_calc_R1, cylinder_calc_err_R1, },
   { "R<sub>2</sub>", N_("Radius 2"), 1, 0, 0, cylinder_calc_R2, cylinder_calc_err_R2, },
};

static const FitShapeParam pring_params[] = {
   { "x<sub>0</sub>", N_("X-coordinate of the center"),     1,  0, 0,                      },
   { "y<sub>0</sub>", N_("Y-coordinate of the center"),     1,  0, 0,                      },
   { "z<sub>0</sub>", N_("Base plane value at the origin"), 0,  1, 0,                      },
   { "R",             N_("Mean circle radius"),             1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "w",             N_("Width at half-height"),           1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "h",             N_("Height of the feature"),          0,  1, 0,                      },
   { "s",             NULL,                                 0,  1, 0,                      },
   { "b<sub>x</sub>", N_("Base plane x-coefficient"),       -1, 1, 0,                      },
   { "b<sub>y</sub>", N_("Base plane y-coefficient"),       -1, 1, 0,                      },
};

#define pring_secondary NULL

static const FitShapeParam lunette_params[] = {
   { "x<sub>0</sub>", N_("X-coordinate of the center"),     1,  0, 0,                      },
   { "y<sub>0</sub>", N_("Y-coordinate of the center"),     1,  0, 0,                      },
   { "z<sub>0</sub>", N_("Base plane value at the origin"), 0,  1, 0,                      },
   { "R",             N_("Outer radius"),                   1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "phi",           N_("Orientation"),                    0,  0, GWY_NLFIT_PARAM_ANGLE,  },
   { "h",             NULL,                                 0,  1, GWY_NLFIT_PARAM_ABSVAL, },
   { "δ",             N_("Distance between centers"),       1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "r",             N_("Inner radius"),                   1,  0, GWY_NLFIT_PARAM_ABSVAL, },
   { "b<sub>x</sub>", N_("Base plane x-coefficient"),       -1, 1, 0,                      },
   { "b<sub>y</sub>", N_("Base plane y-coefficient"),       -1, 1, 0,                      },
   { "p",             NULL,                                 0,  0, GWY_NLFIT_PARAM_ABSVAL, },
};

#define lunette_secondary NULL

static void               enforce_param_flags    (const FitShapeFunc *builtin,
                                                  gdouble *param);
static gdouble            enforce_secondary_flags(const FitShapeFunc *builtin,
                                                  guint i,
                                                  gdouble value);
static GwyShapeFitPreset* create_static_preset   (const FitShapeFunc *function);
static GwySurface*        reduce_data_size       (const GwyXYZ *xyzsrc,
                                                  guint nsrc,
                                                  guint nred);

G_DEFINE_TYPE(GwyShapeFitPreset, gwy_shape_fit_preset, GWY_TYPE_RESOURCE)

static const FitShapeFunc functions[] = {
    { N_("Plane"),                    FALSE, SHAPE_FDIF_ITEM(plane),       },
    { N_("Step (one-sided)"),         FALSE, SHAPE_FDIF_ITEM(step1),       },
    { N_("Step (two-sided)"),         FALSE, SHAPE_FUNC_ITEM(step2),       },
    { N_("Grating (simple)"),         FALSE, SHAPE_FUNC_ITEM(grating),     },
    { N_("Grating (3-level)"),        FALSE, SHAPE_FUNC_ITEM(grating3),    },
    { N_("Holes"),                    FALSE, SHAPE_FUNC_ITEM(holes),       },
    { N_("Sphere"),                   TRUE,  SHAPE_FUNC_ITEM(sphere),      },
    { N_("Half-sphere"),              TRUE,  SHAPE_FUNC_ITEM(hsphere),     },
    { N_("Cylinder (lying)"),         TRUE,  SHAPE_FUNC_ITEM(hcylinder),   },
    { N_("Gaussian"),                 FALSE, SHAPE_FDIF_ITEM(gaussian),    },
    { N_("Lorentzian"),               FALSE, SHAPE_FDIF_ITEM(lorentzian),  },
    { N_("Pyramid (rectangle)"),      FALSE, SHAPE_FUNC_ITEM(pyramid),     },
    { N_("Pyramid (diamond)"),        FALSE, SHAPE_FUNC_ITEM(pyramidx),    },
    { N_("Pyramid (3-sided)"),        FALSE, SHAPE_FUNC_ITEM(pyramid3),    },
    { N_("Parabolic bump"),           FALSE, SHAPE_FUNC_ITEM(parbump),     },
    { N_("Cone"),                     FALSE, SHAPE_FUNC_ITEM(cone),        },
    { N_("Smooth cone"),              TRUE,  SHAPE_FUNC_ITEM(smthcone),    },
    { N_("Smooth pyramid"),           TRUE,  SHAPE_FUNC_ITEM(smthpyr4),    },
    { N_("Smooth pyramid (3-sided)"), TRUE,  SHAPE_FUNC_ITEM(smthpyr3),    },
    { N_("Cylinder"),                 FALSE, SHAPE_FDIF_ITEM(cylinder),    },
    { N_("Ring"),                     FALSE, SHAPE_FUNC_ITEM(pring),       },
    { N_("Lunette"),                  FALSE, SHAPE_FUNC_ITEM(lunette),     },
    { N_("Exponential"),              FALSE, SHAPE_FDIF_ITEM(exponential), },
};

static void
gwy_shape_fit_preset_class_init(GwyShapeFitPresetClass *klass)
{
    GwyResourceClass *parent_class, *res_class = GWY_RESOURCE_CLASS(klass);

    g_type_class_add_private(klass, sizeof(ShapeFitPreset));

    parent_class = GWY_RESOURCE_CLASS(gwy_shape_fit_preset_parent_class);
    res_class->item_type = *gwy_resource_class_get_item_type(parent_class);

    res_class->item_type.type = G_TYPE_FROM_CLASS(klass);
    res_class->name = "shapefitpresets";
    res_class->inventory = gwy_inventory_new(&res_class->item_type);
    gwy_inventory_forget_order(res_class->inventory);
}

static void
gwy_shape_fit_preset_init(GwyShapeFitPreset *preset)
{
    preset->priv = G_TYPE_INSTANCE_GET_PRIVATE(preset, GWY_TYPE_SHAPE_FIT_PRESET, ShapeFitPreset);
}

void
_gwy_shape_fit_preset_class_setup_presets(void)
{
    GwyResourceClass *klass;
    GwyShapeFitPreset *preset;
    guint i;

    /* Force class instantiation, this function is called before it's first referenced. */
    klass = g_type_class_ref(GWY_TYPE_SHAPE_FIT_PRESET);

    for (i = 0; i < G_N_ELEMENTS(functions); i++) {
        preset = create_static_preset(functions + i);
        gwy_inventory_insert_item(klass->inventory, preset);
        g_object_unref(preset);
    }
    gwy_inventory_restore_order(klass->inventory);
    gwy_inventory_set_default_item_name(klass->inventory, functions[0].name);

    /* The presets added a reference so we can safely unref it again */
    g_type_class_unref(klass);
}

static guint
calculate_nreduced(guint n)
{
    return (guint)sqrt(n*(gdouble)NREDLIM);
}

/**
 * gwy_shape_fit_preset_needs_same_units:
 * @preset: A 3D geometrical shape fitting function.
 *
 * Reports if a 3D geometrical shape fitter preset requires the same lateral and value units.
 *
 * For instance, fitting a sphere is meaningless if the horizontal and vertical radii would be different physical
 * quantities.
 *
 * Returns: %TRUE if the function requires the same lateral and value units.
 *
 * Since: 2.47
 **/
gboolean
gwy_shape_fit_preset_needs_same_units(GwyShapeFitPreset *preset)
{
    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), FALSE);
    return preset->priv->builtin->needs_same_units;
}

/**
 * gwy_shape_fit_preset_get_nparams:
 * @preset: A 3D geometrical shape fitting function.
 *
 * Reports the number of parameters of a 3D geometrical shape fitter preset.
 *
 * Returns: The number of function parameters.
 *
 * Since: 2.47
 **/
guint
gwy_shape_fit_preset_get_nparams(GwyShapeFitPreset* preset)
{
    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0);
    return preset->priv->builtin->nparams;
}

/**
 * gwy_shape_fit_preset_get_param_name:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Parameter number.
 *
 * Gets the name of a fitting parameter of a 3D geometrical shape fitter preset.
 *
 * The name may contain Pango markup.
 *
 * Returns: The name of the @i-th parameter.
 *
 * Since: 2.47
 **/
const gchar*
gwy_shape_fit_preset_get_param_name(GwyShapeFitPreset* preset,
                                    guint i)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), "");
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nparams, "");

    return builtin->param[i].name;
}

/**
 * gwy_shape_fit_preset_get_param_description:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Parameter number.
 *
 * Gets the description of a fitting parameter of a 3D geometrical shape fitter preset.
 *
 * The description may contain Pango markup.  It may also be %NULL if the parameter is currently undescribed.
 *
 * Returns: The description of the @i-th parameter, possibly %NULL.
 *
 * Since: 2.61
 **/
const gchar*
gwy_shape_fit_preset_get_param_description(GwyShapeFitPreset *preset,
                                           guint i)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), NULL);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nparams, NULL);

    return builtin->param[i].desc;
}

/**
 * gwy_shape_fit_preset_get_param_flags:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Parameter number.
 *
 * Gets the properties of a fitting parameter of a 3D geometrical shape fitter preset.
 *
 * Returns: The flags of the @i-th parameter.
 *
 * Since: 2.47
 **/
GwyNLFitParamFlags
gwy_shape_fit_preset_get_param_flags(GwyShapeFitPreset* preset,
                                     guint i)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nparams, 0);

    return builtin->param[i].flags;
}

/**
 * gwy_shape_fit_preset_get_param_units:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Parameter number.
 * @siunit_xy: SI unit of lateral coordinates.
 * @siunit_z: SI unit of values.
 *
 * Derives the SI unit of a fitting parameter from the units of abscissa and ordinate.
 *
 * Note that angle parameters are by default in radians and thus unitless. If you want to convert them to degrees for
 * presentation to the user you must do it explicitly.
 *
 * Returns: A newly created #GwySIUnit with the units of the @i-th parameter.
 *
 * Since: 2.47
 **/
GwySIUnit*
gwy_shape_fit_preset_get_param_units(GwyShapeFitPreset *preset,
                                     guint i,
                                     GwySIUnit *siunit_xy,
                                     GwySIUnit *siunit_z)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), NULL);
    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit_xy), NULL);
    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit_z), NULL);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nparams, NULL);

    return gwy_si_unit_power_multiply(siunit_xy, builtin->param[i].power_xy,
                                      siunit_z, builtin->param[i].power_z,
                                      NULL);
}

/**
 * gwy_shape_fit_preset_get_param_power_xy:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Parameter number.
 *
 * Gets the power of abscissa units in a shape fitting parameter.
 *
 * Usually it is easier to let gwy_shape_fit_preset_get_param_units() derive the parameter units.
 *
 * Returns: The power of abscissa in the parameter.
 *
 * Since: 2.50
 **/
gint
gwy_shape_fit_preset_get_param_power_xy(GwyShapeFitPreset *preset,
                                        guint i)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nparams, 0);
    return builtin->param[i].power_xy;
}

/**
 * gwy_shape_fit_preset_get_param_power_z:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Parameter number.
 *
 * Gets the power of ordinate units in a shape fitting parameter.
 *
 * Usually it is easier to let gwy_shape_fit_preset_get_param_units() derive the parameter units.
 *
 * Returns: The power of ordinate in the parameter.
 *
 * Since: 2.50
 **/
gint
gwy_shape_fit_preset_get_param_power_z(GwyShapeFitPreset *preset,
                                       guint i)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nparams, 0);
    return builtin->param[i].power_z;
}

/**
 * gwy_shape_fit_preset_get_nsecondary:
 * @preset: A 3D geometrical shape fitting function.
 *
 * Reports the number of secondary (derived) quantities of a 3D geometrical shape fitter preset.
 *
 * Returns: The number of secondary quantities.
 *
 * Since: 2.47
 **/
guint
gwy_shape_fit_preset_get_nsecondary(GwyShapeFitPreset *preset)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0);
    builtin = preset->priv->builtin;
    return builtin->nsecondary;
}

/**
 * gwy_shape_fit_preset_get_secondary_name:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Secondary quantity number.
 *
 * Gets the name of a secondary (derived) quantity of a 3D geometrical shape fitter preset.
 *
 * The name may contain Pango markup.
 *
 * Returns: The name of the @i-th secondary quantity.
 *
 * Since: 2.47
 **/
const gchar*
gwy_shape_fit_preset_get_secondary_name(GwyShapeFitPreset *preset,
                                        guint i)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), "");
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nsecondary, "");
    return builtin->secondary[i].name;
}

/**
 * gwy_shape_fit_preset_get_secondary_description:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Secondary quantity number.
 *
 * Gets the description of a secondary quantity of a 3D geometrical shape fitter preset.
 *
 * The description may contain Pango markup.  It may also be %NULL if the quantity is currently undescribed.
 *
 * Returns: The description of the @i-th secondary quantity, possibly %NULL.
 *
 * Since: 2.61
 **/
const gchar*
gwy_shape_fit_preset_get_secondary_description(GwyShapeFitPreset *preset,
                                               guint i)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), "");
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nsecondary, "");
    return builtin->secondary[i].desc;
}

/**
 * gwy_shape_fit_preset_get_secondary_flags:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Secondary quantity number.
 *
 * Gets the properties of a secondary (derived) quantity of a 3D geometrical shape fitter preset.
 *
 * Returns: The flags of the @i-th secondary quantity.
 *
 * Since: 2.47
 **/
GwyNLFitParamFlags
gwy_shape_fit_preset_get_secondary_flags(GwyShapeFitPreset *preset,
                                         guint i)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nsecondary, 0);
    return builtin->secondary[i].flags;
}

/**
 * gwy_shape_fit_preset_get_secondary_value:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Secondary quantity number.
 * @param: Values of fitting parameters for function @preset.
 *
 * Calculates the value of a secondary (derived) quantity of a 3D geometrical shape fitter preset.
 *
 * Returns: The value of the @i-th secondary quantity.
 *
 * Since: 2.47
 **/
gdouble
gwy_shape_fit_preset_get_secondary_value(GwyShapeFitPreset *preset,
                                         guint i,
                                         const gdouble *param)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0.0);
    g_return_val_if_fail(param, 0.0);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nsecondary, 0.0);
    return enforce_secondary_flags(builtin, i, builtin->secondary[i].calc(param));
}

/**
 * gwy_shape_fit_preset_get_secondary_error:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Secondary quantity number.
 * @param: Values of fitting parameters for function @preset.
 * @error: Values of errors of fitting parameters for function @preset.
 * @correl: Parameter correlation matrix for function @preset (passed as lower triangular matrix).
 *
 * Calculates the error of a secondary (derived) quantity of a 3D geometrical shape fitter preset.
 *
 * The error is calculated by numerical differentiation of the function and applying the law of error propagation.
 *
 * Returns: The error of the @i-th secondary quantity.
 *
 * Since: 2.47
 **/
gdouble
gwy_shape_fit_preset_get_secondary_error(GwyShapeFitPreset *preset,
                                         guint i,
                                         const gdouble *param,
                                         const gdouble *error,
                                         const gdouble *correl)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0.0);
    g_return_val_if_fail(param, 0.0);
    g_return_val_if_fail(error, 0.0);
    g_return_val_if_fail(correl, 0.0);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nsecondary, 0.0);
    return builtin->secondary[i].calc_err(param, error, correl);
}

/**
 * gwy_shape_fit_preset_get_secondary_units:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Secondary quantity number.
 * @siunit_xy: SI unit of lateral coordinates.
 * @siunit_z: SI unit of values.
 *
 * Derives the SI unit of a secondary (derived) quantity from the units of abscissa and ordinate.
 *
 * Note that angle parameters are by default in radians and thus unitless. If you want to convert them to degrees for
 * presentation to the user you must do it explicitly.
 *
 * Returns: A newly created #GwySIUnit with the units of the @i-th secondary quantity.
 *
 * Since: 2.47
 **/
GwySIUnit*
gwy_shape_fit_preset_get_secondary_units(GwyShapeFitPreset *preset,
                                         guint i,
                                         GwySIUnit *siunit_xy,
                                         GwySIUnit *siunit_z)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), NULL);
    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit_xy), NULL);
    g_return_val_if_fail(GWY_IS_SI_UNIT(siunit_z), NULL);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nsecondary, NULL);

    return gwy_si_unit_power_multiply(siunit_xy, builtin->secondary[i].power_xy,
                                      siunit_z, builtin->secondary[i].power_z,
                                      NULL);
}

/**
 * gwy_shape_fit_preset_get_secondary_power_xy:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Secondary quantity number.
 *
 * Gets the power of abscissa units in a shape fitting secondary (derived) quantity.
 *
 * Usually it is easier to let gwy_shape_fit_preset_get_secondary_units() derive the quantity units.
 *
 * Returns: The power of abscissa in the secondary quantity.
 *
 * Since: 2.50
 **/
gint
gwy_shape_fit_preset_get_secondary_power_xy(GwyShapeFitPreset *preset,
                                            guint i)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nsecondary, 0);
    return builtin->secondary[i].power_xy;
}

/**
 * gwy_shape_fit_preset_get_secondary_power_z:
 * @preset: A 3D geometrical shape fitting function.
 * @i: Secondary quantity number.
 *
 * Gets the power of ordinate units in a shape fitting secondary (derived) quantity.
 *
 * Usually it is easier to let gwy_shape_fit_preset_get_secondary_units() derive the quantity units.
 *
 * Returns: The power of ordinate in the secondary quantity.
 *
 * Since: 2.50
 **/
gint
gwy_shape_fit_preset_get_secondary_power_z(GwyShapeFitPreset *preset,
                                           guint i)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0);
    builtin = preset->priv->builtin;
    g_return_val_if_fail(i < builtin->nsecondary, 0);
    return builtin->secondary[i].power_z;
}

/**
 * gwy_shape_fit_preset_setup:
 * @preset: A 3D geometrical shape fitting function.
 * @points: Array of XYZ data to fit.
 * @n: Number of data points.
 * @params: The array to fill with initialised parameter values.
 *
 * Initialises parameter values of a 3D geometrical shape fitter preset.
 *
 * The parameters are quickly set to reasonable values that roughly correspond to the ranges of the data points.  They
 * may serve as starting values for manual experimentation but often will not be good enough as initial parameter
 * estimates for the fit.  See also gwy_shape_fit_preset_guess().
 *
 * Since: 2.47
 **/
void
gwy_shape_fit_preset_setup(GwyShapeFitPreset *preset,
                           const GwyXYZ *points,
                           guint n,
                           gdouble *params)
{
    FitShapeEstimateCache estimcache;
    const FitShapeFunc *builtin;

    g_return_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset));
    g_return_if_fail(points);
    g_return_if_fail(params);
    builtin = preset->priv->builtin;

    gwy_clear(&estimcache, 1);
    /* The function has a return value but it should always succeed. */
    builtin->initialise(points, n, params, &estimcache);
    enforce_param_flags(builtin, params);
}

/**
 * gwy_shape_fit_preset_guess:
 * @preset: A 3D geometrical shape fitting function.
 * @points: Array of XYZ data to fit.
 * @n: Number of data points.
 * @params: The array to fill with initialised parameter values.
 *
 * Estimates parameter values of a 3D geometrical shape fitter preset.
 *
 * This function tries to find initial parameter estimates that are good enough for the fit the converge.  Of course,
 * it is not guaranteed it always succeeds.  For some shapes it can be noticeably slower than
 * gwy_shape_fit_preset_setup().
 *
 * The estimate may not be deterministic.  For large point sets some estimates are carried out using a randomly
 * selected subset of points.
 *
 * If the function cannot find how the data points could correspond to the preset geometrical shape it return %FALSE.
 * Parameter values are still set. However, in this case they may be no better than from gwy_shape_fit_preset_setup().
 *
 * Returns: %TRUE if the estimation succeeded, %FALSE if it failed.
 *
 * Since: 2.47
 **/
gboolean
gwy_shape_fit_preset_guess(GwyShapeFitPreset *preset,
                           const GwyXYZ *points,
                           guint n,
                           gdouble *params)
{
    FitShapeEstimateCache estimcache;
    const FitShapeFunc *builtin;
    gboolean ok;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), FALSE);
    g_return_val_if_fail(points, FALSE);
    g_return_val_if_fail(params, FALSE);
    builtin = preset->priv->builtin;

    gwy_clear(&estimcache, 1);
    ok = builtin->estimate(points, n, params, &estimcache);
    enforce_param_flags(builtin, params);

    return ok;
}

/**
 * gwy_shape_fit_preset_get_value:
 * @preset: A 3D geometrical shape fitting function.
 * @x: X-coordinate.
 * @y: Y-coordinate.
 * @params: Fitting parameter values.
 *
 * Calculates the value of a 3D geometrical shape fitter preset in a single point.
 *
 * If you want multiple values you should use either gwy_shape_fit_preset_calculate_z() or
 * gwy_shape_fit_preset_calculate_xyz() instead of calling this function in a cycle.
 *
 * Returns: The calculated function value in (@x,@y).
 *
 * Since: 2.47
 **/
gdouble
gwy_shape_fit_preset_get_value(GwyShapeFitPreset *preset,
                               gdouble x,
                               gdouble y,
                               const gdouble *params)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), 0.0);
    g_return_val_if_fail(params, 0.0);
    builtin = preset->priv->builtin;
    return builtin->function(x, y, params);
}

/**
 * gwy_shape_fit_preset_calculate_z:
 * @preset: A 3D geometrical shape fitting function.
 * @points: Array of @n XYZ data defining the lateral coordinates.
 * @z: Array length @n to fill with calculated values.
 * @n: Number of items in @points and @z.
 * @params: Fitting parameter values.
 *
 * Calculates values of a 3D geometrical shape fitter preset in an array of points.
 *
 * The z-coordinates in @points are ignored.  Only the lateral coordinates are used.
 *
 * See also gwy_shape_fit_preset_calculate_xyz().
 *
 * Since: 2.47
 **/
void
gwy_shape_fit_preset_calculate_z(GwyShapeFitPreset *preset,
                                 const GwyXYZ *points,
                                 gdouble *z,
                                 guint n,
                                 const gdouble *params)
{
    const FitShapeFunc *builtin;
    FitShapeXYFunc func;
    guint i;

    g_return_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset));
    g_return_if_fail(params);
    if (!n)
        return;

    g_return_if_fail(points);
    g_return_if_fail(z);
    builtin = preset->priv->builtin;
    func = builtin->function;
    for (i = 0; i < n; i++)
        z[i] = func(points[i].x, points[i].y, params);
}

/**
 * gwy_shape_fit_preset_calculate_xyz:
 * @preset: A 3D geometrical shape fitting function.
 * @points: Array of @n XYZ data defining the lateral coordinates.  The z-coordinates will be filled with the
 *          calculated values.
 * @n: Number of items in @points.
 * @params: Fitting parameter values.
 *
 * Calculates values of a 3D geometrical shape fitter preset in an array of points.
 *
 * See also gwy_shape_fit_preset_calculate_z().
 *
 * Since: 2.47
 **/
void
gwy_shape_fit_preset_calculate_xyz(GwyShapeFitPreset *preset,
                                   GwyXYZ *points,
                                   guint n,
                                   const gdouble *params)
{
    const FitShapeFunc *builtin;
    FitShapeXYFunc func;
    guint i;

    g_return_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset));
    g_return_if_fail(params);
    if (!n)
        return;

    g_return_if_fail(points);
    builtin = preset->priv->builtin;
    func = builtin->function;
    for (i = 0; i < n; i++)
        points[i].z = func(points[i].x, points[i].y, params);
}

/**
 * gwy_shape_fit_preset_create_fitter:
 * @preset: A 3D geometrical shape fitting function.
 *
 * Creates a non-linear least-squares fitter for a 3D geometrical shape.
 *
 * The created fitter will be of the opaque indexed data type, as created with gwy_math_nlfit_new_idx().
 *
 * If you do not need to modify the fitter settings you can use gwy_shape_fit_preset_fit() directly with %NULL fitter.
 *
 * Returns: A newly created fitter for @preset.
 *
 * Since: 2.47
 **/
GwyNLFitter*
gwy_shape_fit_preset_create_fitter(GwyShapeFitPreset *preset)
{
    const FitShapeFunc *builtin;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), NULL);
    builtin = preset->priv->builtin;
    return gwy_math_nlfit_new_idx(builtin->fit_function, builtin->diff_function);
}

/**
 * gwy_shape_fit_preset_fit:
 * @preset: A 3D geometrical shape fitting function.
 * @fitter: A Marquardt-Levenberg nonlinear fitter already initialized for @preset's function, or %NULL.
 * @points: Array of @n XYZ data defining the lateral coordinates and values to fit.
 * @n: Number of items in @points.
 * @params: Fitting parameters filled with initial estimates (the fitting starts from the provided values).
 * @fixed_param: Which parameters should be treated as fixed (set corresponding element to %TRUE for them).  May be
 *               %NULL if all parameters are free.
 * @rss: Location to store the residual sum of squares, as returned by gwy_math_nlfit_fit_idx(), may be %NULL.
 *
 * Performs a non-linear least-squares fit with a 3D geometrical shape fitter.
 *
 * If you pass %NULL @fitter the function creates one for you and immediately performs the fit.  If you want to modify
 * the fitter settings beforehand or set callback functions create it using gwy_shape_fit_preset_create_fitter() and
 * pass to this function.  The fitter must be created for the same preset.
 *
 * Additional quantities such as parameter errors or the correlation matrix can be obtained from the fitter. See
 * gwy_math_nlfit_fit_full() for details.
 *
 * Returns: Either @fitter itself, or a newly created fitter if it was %NULL.
 *
 * Since: 2.47
 **/
GwyNLFitter*
gwy_shape_fit_preset_fit(GwyShapeFitPreset *preset,
                         GwyNLFitter *fitter,
                         const GwyXYZ *points,
                         guint n,
                         gdouble *params,
                         const gboolean *fixed_param,
                         gdouble *rss)
{
    ShapeFitPreset *priv;
    const FitShapeFunc *builtin;
    gdouble myrss;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), NULL);
    g_return_val_if_fail(points, NULL);
    g_return_val_if_fail(params, NULL);
    if (!fitter)
        fitter = gwy_shape_fit_preset_create_fitter(preset);

    priv = preset->priv;
    priv->xyz = points;
    builtin = priv->builtin;
    myrss = gwy_math_nlfit_fit_idx_full(fitter, n, builtin->nparams, params, fixed_param, NULL, priv);
    priv->xyz = NULL;
    enforce_param_flags(builtin, params);
    if (rss)
        *rss = myrss;

    return fitter;
}

/**
 * gwy_shape_fit_preset_quick_fit:
 * @preset: A 3D geometrical shape fitting function.
 * @fitter: A Marquardt-Levenberg nonlinear fitter already initialized for @preset's function, or %NULL.
 * @points: Array of @n XYZ data defining the lateral coordinates and values to fit.
 * @n: Number of items in @points.
 * @params: Fitting parameters filled with initial estimates (the fitting starts from the provided values).
 * @fixed_param: Which parameters should be treated as fixed (set corresponding element to %TRUE for them).  May be
 *               %NULL if all parameters are free.
 * @rss: Location to store the residual sum of squares, as returned by gwy_math_nlfit_fit_idx(), may be %NULL.
 *
 * Performs a rough non-linear least-squares fit with a 3D geometrical shape fitter.
 *
 * See gwy_shape_fit_preset_fit() for discussion.  This functions differs by using a reduced number of points to
 * perform the fit (unless the input set is small) and may also set fitter parameters to make the least squares method
 * terminate faster.  The input points are reduced the same way as in gwy_surface_reduce_points().
 *
 * Returns: Either @fitter itself, or a newly created fitter if it was %NULL.
 *
 * Since: 2.48
 **/
GwyNLFitter*
gwy_shape_fit_preset_quick_fit(GwyShapeFitPreset *preset,
                               GwyNLFitter *fitter,
                               const GwyXYZ *points,
                               guint n,
                               gdouble *params,
                               const gboolean *fixed_param,
                               gdouble *rss)
{
    ShapeFitPreset *priv;
    const FitShapeFunc *builtin;
    guint nred = calculate_nreduced(n);
    gdouble myrss;

    g_return_val_if_fail(GWY_IS_SHAPE_FIT_PRESET(preset), NULL);
    g_return_val_if_fail(points, NULL);
    g_return_val_if_fail(params, NULL);
    if (!fitter)
        fitter = gwy_shape_fit_preset_create_fitter(preset);

    gwy_math_nlfit_set_max_iterations(fitter, 30);
    priv = preset->priv;
    builtin = priv->builtin;
    if (nred < n) {
        GwySurface *surface = reduce_data_size(points, n, nred);
        priv->xyz = surface->data;
        myrss = gwy_math_nlfit_fit_idx_full(fitter, nred, builtin->nparams, params, fixed_param, NULL, priv);
        g_object_unref(surface);
    }
    else {
        priv->xyz = points;
        myrss = gwy_math_nlfit_fit_idx_full(fitter, n, builtin->nparams, params, fixed_param, NULL, priv);
    }
    priv->xyz = NULL;
    enforce_param_flags(builtin, params);
    if (rss)
        *rss = myrss;

    return fitter;
}

static void
enforce_param_flags(const FitShapeFunc *builtin, gdouble *param)
{
    guint i;

    for (i = 0; i < builtin->nparams; i++) {
        GwyNLFitParamFlags flags = builtin->param[i].flags;

        if (flags & GWY_NLFIT_PARAM_ANGLE)
            param[i] = fmod(param[i], 2.0*G_PI);
        if (flags & GWY_NLFIT_PARAM_ABSVAL)
            param[i] = fabs(param[i]);
    }
}

static gdouble
enforce_secondary_flags(const FitShapeFunc *builtin, guint i, gdouble value)
{
    GwyNLFitParamFlags flags;

    g_return_val_if_fail(i < builtin->nsecondary, value);
    flags = builtin->secondary[i].flags;

    if (flags & GWY_NLFIT_PARAM_ANGLE)
        value = fmod(value, 2.0*G_PI);
    if (flags & GWY_NLFIT_PARAM_ABSVAL)
        value = fabs(value);

    return value;
}

static GwySurface*
reduce_data_size(const GwyXYZ *xyzsrc, guint nsrc, guint nred)
{
    GwySurface *reduced, *surface = gwy_surface_new();
    guint n = surface->n;
    GwyXYZ *data = surface->data;

    /* XXX: This is a dirty trick we can only dare to attempt since we are doing it inside the library and know
     * gwy_surface_reduce_points() does not change the @data and @n fields.  Do not try it at home! */
    surface->data = (GwyXYZ*)xyzsrc;
    surface->n = nsrc;
    reduced = gwy_surface_reduce_points(surface, nred);
    surface->data = data;
    surface->n = n;

    return reduced;
}

static GwyShapeFitPreset*
create_static_preset(const FitShapeFunc *builtin)
{
    GwyShapeFitPreset *preset;

    preset = g_object_new(GWY_TYPE_SHAPE_FIT_PRESET, "is-const", TRUE, NULL);
    preset->priv->builtin = builtin;
    g_string_assign(GWY_RESOURCE(preset)->name, builtin->name);

    return preset;
}

/**
 * gwy_shape_fit_presets:
 *
 * Gets inventory with all the 3D geometric shape fitting presets.
 *
 * Returns: 3D geometric shape fitting preset inventory.
 *
 * Since: 2.47
 **/
GwyInventory*
gwy_shape_fit_presets(void)
{
    GTypeClass *klass = g_type_class_peek(GWY_TYPE_SHAPE_FIT_PRESET);
    return GWY_RESOURCE_CLASS(klass)->inventory;
}

/*******************************************************************************************************************
 *
 * General estimator helpers and math support functions.
 *
 *******************************************************************************************************************/

#ifdef HAVE_SINCOS
#define _gwy_sincos sincos
#else
static inline void
_gwy_sincos(gdouble x, gdouble *s, gdouble *c)
{
    *s = sin(x);
    *c = cos(x);
}
#endif

/* cosh(x) - 1, safe for small arguments */
static inline gdouble
gwy_coshm1(gdouble x)
{
    gdouble x2 = x*x;
    if (x2 > 3e-5)
        return cosh(x) - 1.0;
    return x2*(0.5 + x2/24.0);
}

static inline gdouble
diff_numerically(FitShapeXYFunc func, gdouble x, gdouble y, gdouble f0,
                 gdouble *param, guint i,
                 gdouble nat_dp, gboolean use_nat_dp)
{
    gdouble p0, dp, der;

    if (use_nat_dp || !(param[i] > 0.0))
        dp = nat_dp;
    else
        dp = 1e-5*fabs(param[i]);

    p0 = param[i];
    param[i] += dp;
    der = (func(x, y, param) - f0)/dp;
    param[i] = p0;

    return der;
}

#if 0
static void
compare_der(gdouble x, gdouble y,
            const gdouble *param, const gdouble *der, guint nparam,
            FitShapeXYFunc func, guint is_good)
{
    gdouble f0, p[nparam];
    guint i;

    g_printerr("%g %g", x, y);
    gwy_assign(p, param, nparam);
    f0 = func(x, y, p);
    for (i = 0; i < nparam; i++) {
        if (is_good & (1 << i))
            continue;

        g_printerr(" (%u: %g %g)",
                   i, der[i],
                   diff_numerically(func, x, y, f0, p, i, 1e-5, FALSE));
    }
    g_printerr("\n");
}
#endif

/* Mean value of xy point cloud (not necessarily centre, that depends on
 * the density). */
static void
mean_x_y(const GwyXYZ *xyz, guint n, gdouble *pxm, gdouble *pym,
         FitShapeEstimateCache *estimcache)
{
    gdouble xm = 0.0, ym = 0.0;
    guint i;

    if (estimcache && estimcache->have_mean) {
        gwy_debug("using cache %p", estimcache);
        *pxm = estimcache->xm;
        *pym = estimcache->ym;
        return;
    }

    if (!n) {
        *pxm = *pym = 0.0;
        return;
    }

    for (i = 0; i < n; i++) {
        xm += xyz[i].x;
        ym += xyz[i].y;
    }

    *pxm = xm/n;
    *pym = ym/n;

    if (estimcache) {
        gwy_debug("filling cache %p", estimcache);
        estimcache->have_mean = TRUE;
        estimcache->xm = *pxm;
        estimcache->ym = *pym;
    }
}

/* Minimum and maximum of an array of values. */
static void
range_z(const GwyXYZ *xyz, guint n, gdouble *pmin, gdouble *pmax,
        FitShapeEstimateCache *estimcache)
{
    gdouble min, max;
    guint i;

    if (estimcache && estimcache->have_zrange) {
        gwy_debug("using cache %p", estimcache);
        *pmin = estimcache->zmin;
        *pmax = estimcache->zmax;
        return;
    }

    if (!n) {
        *pmin = *pmax = 0.0;
        return;
    }

    min = max = xyz[0].z;
    for (i = 1; i < n; i++) {
        if (xyz[i].z < min)
            min = xyz[i].z;
        if (xyz[i].z > max)
            max = xyz[i].z;
    }

    *pmin = min;
    *pmax = max;

    if (estimcache) {
        gwy_debug("filling cache %p", estimcache);
        estimcache->have_zrange = TRUE;
        estimcache->zmin = *pmin;
        estimcache->zmax = *pmax;
    }
}

/* Simple stats of an array of values. */
static void
stat_z(const GwyXYZ *xyz, guint n,
       gdouble *zmean, gdouble *zrms, gdouble *zskew,
       FitShapeEstimateCache *estimcache)
{
    gdouble s = 0.0, s2 = 0.0, s3 = 0.0;
    guint i;

    if (estimcache && estimcache->have_zstats) {
        gwy_debug("using cache %p", estimcache);
        if (zmean)
            *zmean = estimcache->zmean;
        if (zrms)
            *zrms = estimcache->zrms;
        if (zskew)
            *zskew = estimcache->zskew;
        return;
    }

    if (!n) {
        if (zmean)
            *zmean = 0.0;
        if (zrms)
            *zrms = 0.0;
        if (zskew)
            *zskew = 0.0;
        return;
    }

    for (i = 0; i < n; i++)
        s += xyz[i].z;
    s /= n;

    for (i = 0; i < n; i++) {
        gdouble d = xyz[i].z - s;
        s2 += d*d;
        s3 += d*d*d;
    }

    if (s2) {
        s2 = sqrt(s2/n);
        s3 /= n*s2*s2*s2;
    }

    if (zmean)
        *zmean = s;
    if (zrms)
        *zrms = s2;
    if (zskew)
        *zskew = s3;

    if (estimcache) {
        gwy_debug("filling cache %p", estimcache);
        estimcache->have_zstats = TRUE;
        estimcache->zmean = s;
        estimcache->zrms = s2;
        estimcache->zskew = s3;
    }
}

/* Approximately cicrumscribe a set of points by finding a containing
 * octagon. */
static void
circumscribe_x_y(const GwyXYZ *xyz, guint n,
                 gdouble *pxc, gdouble *pyc, gdouble *pr,
                 FitShapeEstimateCache *estimcache)
{
    gdouble min[4], max[4], r[4];
    guint i, j;

    if (estimcache && estimcache->have_circle) {
        gwy_debug("using cache %p", estimcache);
        *pxc = estimcache->xc;
        *pyc = estimcache->yc;
        *pr = estimcache->r;
        return;
    }

    if (!n) {
        *pxc = *pyc = 0.0;
        *pr = 1.0;
        return;
    }

    for (j = 0; j < 4; j++) {
        min[j] = G_MAXDOUBLE;
        max[j] = -G_MAXDOUBLE;
    }

    for (i = 0; i < n; i++) {
        gdouble x = xyz[i].x, y = xyz[i].y;
        gdouble t[4] = { x, x+y, y, y-x };

        for (j = 0; j < 4; j++) {
            if (t[j] < min[j])
                min[j] = t[j];
            if (t[j] > max[j])
                max[j] = t[j];
        }
    }

    for (j = 0; j < 4; j++) {
        r[j] = sqrt(10.0)/3.0*(max[j] - min[j]);
        if (j % 2)
            r[j] /= G_SQRT2;
    }

    i = 0;
    for (j = 1; j < 4; j++) {
        if (r[j] > r[i])
            i = j;
    }

    *pr = 0.5*r[i];
    if (i % 2) {
        *pxc = (min[1] - min[3] + max[1] - max[3])/4.0;
        *pyc = (min[1] + min[3] + max[1] + max[3])/4.0;
    }
    else {
        *pxc = (min[0] + max[0])/2.0;
        *pyc = (min[2] + max[2])/2.0;
    }

    if (estimcache) {
        gwy_debug("filling cache %p", estimcache);
        estimcache->have_circle = TRUE;
        estimcache->xc = *pxc;
        estimcache->yc = *pyc;
        estimcache->r = *pr;
    }
}

/* Project xyz point cloud to a line rotated by angle phi anti-clockwise
 * from the horizontal line (x axis). */
static gdouble
projection_to_line(const GwyXYZ *xyz,
                   guint n,
                   gdouble phi,
                   gdouble xc, gdouble yc,
                   GwyDataLine *mean_line,
                   GwyDataLine *rms_line,
                   guint *counts)
{
    guint res = gwy_data_line_get_res(mean_line);
    gdouble *mean = gwy_data_line_get_data(mean_line);
    gdouble *rms = rms_line ? gwy_data_line_get_data(rms_line) : NULL;
    gdouble dx = gwy_data_line_get_real(mean_line)/res;
    gdouble off = gwy_data_line_get_offset(mean_line);
    gdouble c = cos(phi), s = sin(phi), total_ms = 0.0;
    guint i, total_n = 0;
    gint j;

    gwy_data_line_clear(mean_line);
    gwy_clear(counts, res);

    for (i = 0; i < n; i++) {
        gdouble x = xyz[i].x - xc, y = xyz[i].y - yc;
        x = x*c - y*s;
        j = (gint)floor((x - off)/dx);
        if (j >= 0 && j < res) {
            mean[j] += xyz[i].z;
            counts[j]++;
        }
    }

    for (j = 0; j < res; j++) {
        if (counts[j]) {
            mean[j] /= counts[j];
        }
    }

    if (!rms_line)
        return 0.0;

    gwy_data_line_clear(rms_line);

    for (i = 0; i < n; i++) {
        gdouble x = xyz[i].x - xc, y = xyz[i].y - yc;
        x = x*c - y*s;
        j = (gint)floor((x - off)/dx);
        if (j >= 0 && j < res)
            rms[j] += (xyz[i].z - mean[j])*(xyz[i].z - mean[j]);
    }

    for (j = 0; j < res; j++) {
        if (counts[j]) {
            total_ms += rms[j];
            rms[j] = sqrt(rms[j]/counts[j]);
            total_n += counts[j];
        }
    }

    return sqrt(total_ms/total_n);
}

/* Find direction along which projections capture best the shape, i.e. most variance remains in the line-averaged
 * data.  The returned angle is rotation of the axis anti-clockwise with respect to the x-axis. */
static gdouble
estimate_projection_direction(const GwyXYZ *xyz, guint n,
                              FitShapeEstimateCache *estimcache)
{
    enum { NROUGH = 60, NFINE = 8 };

    GwyDataLine *mean_line, *rms_line;
    guint *counts;
    gdouble xc, yc, r, phi, alpha0, alpha_step, rms;
    gdouble best_rms = G_MAXDOUBLE, best_alpha = 0.0;
    guint iter, i, ni, res;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    res = (guint)floor(0.8*sqrt(n) + 1.0);

    mean_line = gwy_data_line_new(res, 2.0*r, FALSE);
    gwy_data_line_set_offset(mean_line, -r);
    rms_line = gwy_data_line_new_alike(mean_line, FALSE);
    counts = g_new(guint, res);

    for (iter = 0; iter < 6; iter++) {
        if (iter == 0) {
            ni = NROUGH;
            alpha_step = G_PI/ni;
            alpha0 = 0.0;
        }
        else {
            /* Choose the fine points so that we do not repeat calculation in any of the rough points. */
            ni = NFINE;
            alpha0 = best_alpha - alpha_step*(NFINE - 1.0)/(NFINE + 1.0);
            alpha_step = 2.0*alpha_step/(NFINE + 1.0);
        }

        for (i = 0; i < ni; i++) {
            phi = alpha0 + i*alpha_step;
            rms = projection_to_line(xyz, n, phi, xc, yc, mean_line, rms_line, counts);
            gwy_debug("[%u] %g %g", iter, phi, rms);
            if (rms < best_rms) {
                best_rms = rms;
                best_alpha = phi;
            }
        }
    }

    g_object_unref(mean_line);
    g_object_unref(rms_line);
    g_free(counts);

    if (best_alpha > 0.5*G_PI)
        best_alpha += G_PI;

    return best_alpha;
}

/* Estimate projection direction, possibly on reduced data.  This is useful when the estimator does not need reduced
 * data for anything else. */
static gdouble
estimate_projection_direction_red(const GwyXYZ *xyz, guint n,
                                  FitShapeEstimateCache *estimcache)
{
    FitShapeEstimateCache estimcachered;
    guint nred = calculate_nreduced(n);
    GwySurface *surface;
    gdouble phi;

    if (nred >= n)
        return estimate_projection_direction(xyz, n, estimcache);

    /* Make sure caching still works for the reduced data. */
    surface = reduce_data_size(xyz, n, nred);
    gwy_clear(&estimcachered, 1);
    phi = estimate_projection_direction(gwy_surface_get_data_const(surface), nred, &estimcachered);
    g_object_unref(surface);

    return phi;
}

static void
data_line_shorten(GwyDataLine *dline, const guint *counts, guint threshold)
{
    guint res = gwy_data_line_get_res(dline);
    guint from = 0, to = res-1;
    gdouble off;

    while (to > from && counts[to] < threshold)
        to--;
    while (from < to && counts[from] < threshold)
        from++;

    off = (from*gwy_data_line_get_real(dline)/res
           + gwy_data_line_get_offset(dline));

    gwy_data_line_resize(dline, from, to+1);
    gwy_data_line_set_offset(dline, off);
}

/* Estimate the period of a periodic structure, knowing already the rotation. The returned phase is such that if you
 * subtract it from the rotated abscissa value then the projection will have a positive peak (some kind of maximum)
 * centered around zero, whatever that means for specific grating-like structures.  */
static gboolean
estimate_period_and_phase(const GwyXYZ *xyz, guint n,
                          gdouble phi, gdouble *pT, gdouble *poff,
                          FitShapeEstimateCache *estimcache)
{
    GwyDataLine *mean_line, *tmp_line;
    gdouble xc, yc, r, T, t, real, off, a_s, a_c, phi0, av, bv;
    const gdouble *mean, *tmp;
    guint *counts;
    guint res, i, ibest;
    gboolean found;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    /* Using more sqrt(n) than can make the sampling too sparse, causing noise and oscillations. */
    res = (guint)floor(0.8*sqrt(n) + 1.0);

    *pT = r/4.0;
    *poff = 0.0;

    mean_line = gwy_data_line_new(res, 2.0*r, FALSE);
    gwy_data_line_set_offset(mean_line, -r);
    tmp_line = gwy_data_line_new_alike(mean_line, FALSE);
    counts = g_new(guint, res);

    projection_to_line(xyz, n, phi, xc, yc, mean_line, NULL, counts);
    data_line_shorten(mean_line, counts, 4);
    g_free(counts);

    res = gwy_data_line_get_res(mean_line);
    gwy_data_line_get_line_coeffs(mean_line, &av, &bv);
    gwy_data_line_line_level(mean_line, av, bv);
    gwy_data_line_psdf(mean_line, tmp_line, GWY_WINDOWING_HANN, GWY_INTERPOLATION_LINEAR);
    tmp = gwy_data_line_get_data_const(tmp_line);

    found = FALSE;
    ibest = G_MAXUINT;
    for (i = 3; i < MIN(res/3, res-3); i++) {
        if (tmp[i] > tmp[i-2] && tmp[i] > tmp[i-1]
            && tmp[i] > tmp[i+1] && tmp[i] > tmp[i+2]) {
            if (ibest == G_MAXUINT || tmp[i] > tmp[ibest]) {
                found = TRUE;
                ibest = i;
            }
        }
    }
    if (!found)
        goto fail;

    T = *pT = 2.0*G_PI/gwy_data_line_itor(tmp_line, ibest);
    gwy_debug("found period %g", T);

    mean = gwy_data_line_get_data_const(mean_line);
    real = gwy_data_line_get_real(mean_line);
    off = gwy_data_line_get_offset(mean_line);
    a_s = a_c = 0.0;
    for (i = 0; i < res; i++) {
        t = off + real/res*(i + 0.5);
        a_s += sin(2*G_PI*t/T)*mean[i];
        a_c += cos(2*G_PI*t/T)*mean[i];
    }
    gwy_debug("a_s %g, a_c %g", a_s, a_c);

    phi0 = atan2(a_s, a_c);
    *poff = phi0*T/(2.0*G_PI) + xc*cos(phi) - yc*sin(phi);

fail:
    g_object_unref(mean_line);
    g_object_unref(tmp_line);

    return found;
}

/* For a shape that consists of a more or less flat base with some feature on it, estimate the base plane (z0) and
 * feature height (h).  The height can be either positive or negative. */
static gboolean
estimate_feature_height(const GwyXYZ *xyz, guint n,
                        gdouble *pz0, gdouble *ph, gdouble *px, gdouble *py,
                        FitShapeEstimateCache *estimcache)
{
    gdouble xm, ym, xc, yc, r, zmin, zmax;
    gdouble r2_large, r2_small;
    gdouble t, zbest, zmean_large = 0.0, zmean_small = 0.0;
    guint i, n_large = 0, n_small = 0;
    gboolean positive;

    if (!n) {
        *pz0 = *ph = 0.0;
        return FALSE;
    }

    range_z(xyz, n, &zmin, &zmax, estimcache);
    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    r2_large = 0.7*r*r;
    r2_small = 0.1*r*r;

    for (i = 0; i < n; i++) {
        gdouble x = xyz[i].x - xc, y = xyz[i].y - yc;
        gdouble r2 = x*x + y*y;

        if (r2 <= r2_small) {
            zmean_small += xyz[i].z;
            n_small++;
        }
        else if (r2 >= r2_large) {
            zmean_large += xyz[i].z;
            n_large++;
        }
    }

    g_assert(n_large);   /* circumscribe_x_y() should ensure this. */
    zmean_large /= n_large;

    if (n_small) {
        zmean_small /= n_small;
        positive = (zmean_small >= zmean_large);
    }
    else
        positive = (fabs(zmean_large - zmin) <= fabs(zmean_large - zmax));

    t = zmax - zmin;
    if (positive) {
        *pz0 = zmin + 0.05*t;
        *ph = 0.9*t;
    }
    else {
        *pz0 = zmax - 0.05*t;
        *ph = -0.9*t;
    }

    xm = 0.0;
    ym = 0.0;
    if (n_small) {
        if (positive) {
            zbest = -G_MAXDOUBLE;
            for (i = 0; i < n; i++) {
                gdouble x = xyz[i].x - xc, y = xyz[i].y - yc;
                gdouble r2 = x*x + y*y;

                if (r2 <= r2_small && xyz[i].z > zbest) {
                    zbest = xyz[i].z;
                    xm = x;
                    ym = y;
                }
            }
        }
        else {
            zbest = G_MAXDOUBLE;
            for (i = 0; i < n; i++) {
                gdouble x = xyz[i].x - xc, y = xyz[i].y - yc;
                gdouble r2 = x*x + y*y;

                if (r2 <= r2_small && xyz[i].z < zbest) {
                    zbest = xyz[i].z;
                    xm = x;
                    ym = y;
                }
            }
        }
    }
    *px = xc + xm;
    *py = yc + ym;

    return TRUE;
}

static gboolean
common_bump_feature_init(const GwyXYZ *xyz, guint n,
                         gdouble *xc, gdouble *yc, gdouble *z0,
                         gdouble *height, gdouble *size,
                         gdouble *a, gdouble *phi,
                         FitShapeEstimateCache *estimcache)
{
    gdouble xm, ym, r, zmin, zmax;

    circumscribe_x_y(xyz, n, &xm, &ym, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    *xc = xm;
    *yc = ym;
    *z0 = zmin;
    *height = zmax - zmin;
    *size = r/3.0;
    *a = 1.0;
    *phi = 0.0;

    return TRUE;
}

static gboolean
common_bump_feature_estimate(const GwyXYZ *xyz, guint n,
                             gdouble *xc, gdouble *yc, gdouble *z0,
                             gdouble *height, gdouble *size,
                             gdouble *a, gdouble *phi,
                             FitShapeEstimateCache *estimcache)
{
    gdouble xm, ym, r;
    guint i, nb;

    /* Just initialise the shape parameters with some sane defaults. */
    *a = 1.0;
    *phi = 0.0;
    circumscribe_x_y(xyz, n, &xm, &ym, &r, estimcache);
    *size = r/3.0;

    if (!estimate_feature_height(xyz, n, z0, height, xc, yc, estimcache))
        return FALSE;
    if (*height == 0.0)
        return TRUE;

    xm = ym = 0.0;
    nb = 0;
    for (i = 0; i < n; i++) {
        gdouble x = xyz[i].x, y = xyz[i].y, z = xyz[i].z;

        if ((z - *z0)/(*height) > 0.4) {
            xm += x;
            ym += y;
            nb++;
        }
    }
    if (!nb)
        return TRUE;

    *xc = xm = xm/nb;
    *yc = ym = ym/nb;

    r = 0.0;
    for (i = 0; i < n; i++) {
        gdouble x = xyz[i].x, y = xyz[i].y, z = xyz[i].z;

        if ((z - *z0)/(*height) > 0.4)
            r += (x - xm)*(x - xm) + (y - ym)*(y - ym);
    }
    *size = sqrt(2.0*r/nb);

    return TRUE;
}

static gdouble
dotprod_with_correl(const gdouble *diff,
                    const gdouble *param_err,
                    const gdouble *correl,
                    guint n)
{
    guint i, j;
    gdouble s = 0.0;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (diff[i] != 0 && diff[j] != 0) {
                gdouble c_ij = (j <= i) ? SLi(correl, i, j) : SLi(correl, j, i);
                gdouble s_ij = c_ij*param_err[i]*param_err[j];
                s += s_ij*diff[i]*diff[j];
            }
        }
    }

    return sqrt(fmax(s, 0.0));
}

static gdouble
data_line_pearson_coeff(GwyDataLine *dline1, GwyDataLine *dline2)
{
    gdouble avg1 = gwy_data_line_get_avg(dline1);
    gdouble avg2 = gwy_data_line_get_avg(dline2);
    gdouble rms1 = gwy_data_line_get_rms(dline1);
    gdouble rms2 = gwy_data_line_get_rms(dline2);
    const gdouble *d1, *d2;
    gdouble c = 0.0;
    guint res, i;

    if (!rms1 || !rms2)
        return 0.0;

    res = gwy_data_line_get_res(dline1);
    g_return_val_if_fail(gwy_data_line_get_res(dline2) == res, 0.0);
    d1 = gwy_data_line_get_data_const(dline1);
    d2 = gwy_data_line_get_data_const(dline2);
    for (i = 0; i < res; i++)
        c += (d1[i] - avg1)*(d2[i] - avg2);

    c /= res*rms1*rms2;
    gwy_debug("%g", c);
    return c;
}

/*******************************************************************************************************************
 *
 * Plane
 *
 *******************************************************************************************************************/

static gdouble
plane_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble z0 = param[0];
    gdouble bx = param[1];
    gdouble by = param[2];

    return z0 + bx*x + by*y;
}

static void
plane_diff(gdouble x, gdouble y,
           G_GNUC_UNUSED const gdouble *param,
           G_GNUC_UNUSED const gboolean *fixed_param,
           gdouble *der)
{
    der[0] = 1;
    der[1] = x;
    der[2] = y;
}

static gboolean
plane_init(const GwyXYZ *xyz, guint n, gdouble *param,
           FitShapeEstimateCache *estimcache)
{
    gdouble zmean;

    stat_z(xyz, n, &zmean, NULL, NULL, estimcache);

    param[0] = zmean;
    param[1] = param[2] = 0.0;

    return TRUE;
}

static gboolean
plane_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
               FitShapeEstimateCache *estimcache)
{
    gdouble xc, yc;
    /* Linear fit with functions 1, x and y. */
    gdouble a[6], b[3];
    guint i;

    /* Using centered coodinates improves the condition number. */
    mean_x_y(xyz, n, &xc, &yc, estimcache);
    gwy_clear(a, 6);
    gwy_clear(b, 3);
    for (i = 0; i < n; i++) {
        gdouble x = xyz[i].x - xc, y = xyz[i].y - yc;

        b[0] += xyz[i].z;
        b[1] += x*xyz[i].z;
        b[2] += y*xyz[i].z;

        a[2] += x*x;
        a[4] += x*y;
        a[5] += y*y;
    }
    a[0] = n;
    a[1] = a[3] = 0.0;

    param[0] = b[0]/n;
    param[1] = 0.0;
    param[2] = 0.0;

    if (!gwy_math_choleski_decompose(3, a))
        return FALSE;

    gwy_math_choleski_solve(3, a, b);

    param[0] = b[0] - xc*b[1] - yc*b[2];
    param[1] = b[1];
    param[2] = b[2];

    return TRUE;
}

static gdouble
plane_calc_phi(const gdouble *param)
{
    return atan2(param[2], param[1]);
}

static gdouble
plane_calc_err_phi(const gdouble *param,
                   const gdouble *param_err,
                   const gdouble *correl)
{
    gdouble diff[G_N_ELEMENTS(sphere_params)];
    gdouble r2 = param[1]*param[1] + param[2]*param[2];

    if (r2 == 0.0)
        return 0.0;

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[1] = -param[2]/r2;
    diff[2] = param[1]/r2;
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
plane_calc_theta(const gdouble *param)
{
    return atan(hypot(param[1], param[2]));
}

static gdouble
plane_calc_err_theta(const gdouble *param,
                     const gdouble *param_err,
                     const gdouble *correl)
{
    gdouble diff[G_N_ELEMENTS(sphere_params)];
    gdouble r2 = param[1]*param[1] + param[2]*param[2];

    if (r2 == 0.0)
        return 0.0;

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[1] = param[1]/(1.0 + r2)/sqrt(r2);
    diff[2] = param[2]/(1.0 + r2)/sqrt(r2);
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

/*******************************************************************************************************************
 *
 * Step (one-sided)
 *
 *******************************************************************************************************************/

static gdouble
step1_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble h = param[0];
    gdouble delta = fabs(param[1]);
    gdouble z0 = param[2];
    gdouble x0 = param[3];
    gdouble phi = param[4];
    gdouble bx = param[5];
    gdouble by = param[6];
    gdouble t, val, cphi, sphi;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi - x0;
    if (G_LIKELY(delta > 0.0))
        val = 0.5*h*(erf(t/delta) + 1.0);
    else {
        if (t < 0.0)
            val = 0.0;
        else if (t > 0.0)
            val = h;
        else
            val = 0.5*h;
    }

    return val + z0 + bx*x + by*y;
}

static void
step1_diff(gdouble x, gdouble y,
           const gdouble *param,
           const gboolean *fixed_param,
           gdouble *der)
{
    gdouble h = param[0];
    gdouble delta = fabs(param[1]);
    /*gdouble z0 = param[2];*/
    gdouble x0 = param[3];
    gdouble phi = param[4];
    /*gdouble bx = param[5];*/
    /*gdouble by = param[6];*/
    gdouble t, cphi, sphi, E;

    der[2] = 1;
    der[5] = x;
    der[6] = y;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi - x0;
    if (G_LIKELY(delta > 0.0)) {
        E = exp(-(t/delta)*(t/delta));

        if (!fixed_param || !fixed_param[0])
            der[0] = 0.5*erf(t/delta) + 0.5;

        der[1] = -h*E/(GWY_SQRT_PI*delta) * t/delta;
        der[3] = -h*E/(GWY_SQRT_PI*delta);
        der[4] = -h*E/(GWY_SQRT_PI*delta) * (x*sphi + y*cphi);
        return;
    }

    /* Derivative by delta is infinity.  Derivatives by x0 and phi are either zero or infinity.  The latter two can
     * give something meaningful on average when done numerically; the derivative by delta cannot. */
    if (t < 0.0) {
        der[0] = 0.0;
        der[1] = 1.0;
    }
    else if (t > 0.0) {
        der[0] = 1.0;
        der[1] = -1.0;
    }
    else {
        der[0] = 0.5;
        der[1] = 0.0;
    }

    {
        gdouble p[G_N_ELEMENTS(step1_params)], f0;

        gwy_assign(p, param, G_N_ELEMENTS(step1_params));
        f0 = step1_func(x, y, p);
        der[3] = diff_numerically(step1_func, x, y, f0, p, 3, 1e-8, FALSE);
        /* Angle has a natural scale. */
        der[4] = diff_numerically(step1_func, x, y, f0, p, 4, 1e-3, TRUE);
    }
}

static gboolean
step1_init(const GwyXYZ *xyz, guint n, gdouble *param,
           FitShapeEstimateCache *estimcache)
{
    gdouble zmin, zmax, xc, yc, r;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    param[0] = zmax - zmin;
    param[1] = 0.01*r;
    param[2] = zmin;
    param[3] = xc;
    param[4] = 0.0;
    param[5] = param[6] = 0.0;

    return TRUE;
}

static gboolean
step1_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
               FitShapeEstimateCache *estimcache)
{
    GwyDataLine *mean_line;
    gdouble zmin, zmax, xc, yc, r, phi, m, dx, off;
    guint *counts = NULL;
    guint i, mpos, res;
    const gdouble *d;
    gboolean inverted;

    res = (guint)floor(0.8*sqrt(n) + 1.0);
    if (res < 8)
        return step1_init(xyz, n, param, estimcache);

    /* First we estimate the orientation (phi). */
    phi = estimate_projection_direction_red(xyz, n, estimcache);
    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    mean_line = gwy_data_line_new(res, 2.0*r, FALSE);
    counts = g_new(guint, res);
    gwy_data_line_set_offset(mean_line, -r);
    projection_to_line(xyz, n, phi, xc, yc, mean_line, NULL, counts);
    data_line_shorten(mean_line, counts, 4);
    d = gwy_data_line_get_data(mean_line);
    dx = gwy_data_line_get_dx(mean_line);
    off = gwy_data_line_get_offset(mean_line);

    mpos = res/2;
    m = 0.0;
    for (i = 1; i < res; i++) {
        if (!counts[i] || !counts[i-1])
           continue;
        if (fabs(d[i] - d[i-1]) > m) {
            m = fabs(d[i] - d[i-1]);
            mpos = i;
        }
    }

    inverted = (d[mpos] - d[mpos-1] < 0.0);
    param[0] = zmax - zmin;
    param[1] = 0.01*r;
    param[2] = zmin;
    param[3] = (inverted ? -1 : 1)*(xc*cos(phi) - yc*sin(phi) + mpos*dx + off);
    param[4] = gwy_canonicalize_angle(phi + (inverted ? G_PI : 0.0), FALSE, TRUE);
    param[5] = param[6] = 0.0;

    return TRUE;
}

/*******************************************************************************************************************
 *
 * Step (two-sided)
 *
 *******************************************************************************************************************/

static gdouble
step2_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble h = param[0];
    gdouble w = fabs(param[1]);
    gdouble delta = fabs(param[2]);
    gdouble z0 = param[3];
    gdouble x0 = param[4];
    gdouble phi = param[5];
    gdouble bx = param[6];
    gdouble by = param[7];
    gdouble t, val, cphi, sphi;

    if (w == 0.0)
        return z0 + bx*x + by*y;

    _gwy_sincos(phi, &sphi, &cphi);
    t = fabs(x*cphi - y*sphi - x0);
    if (G_LIKELY(delta > 0.0))
        val = 0.5*h*(1.0 - erf((t - 0.5*w)/delta));
    else {
        if (t < 0.5*w)
            val = h;
        else if (t > 0.5*w)
            val = 0.0;
        else
            val = 0.5*h;
    }

    return val + z0 + bx*x + by*y;
}

static gboolean
step2_init(const GwyXYZ *xyz, guint n, gdouble *param,
           FitShapeEstimateCache *estimcache)
{
    gdouble zmin, zmax, xc, yc, r;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    param[0] = zmax - zmin;
    param[1] = 0.25*r;
    param[2] = 0.01*r;
    param[3] = zmin;
    param[4] = xc;
    param[5] = 0.0;
    param[6] = param[7] = 0.0;

    return TRUE;
}

static gboolean
step2_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
               FitShapeEstimateCache *estimcache)
{
    GwyDataLine *mean_line;
    gdouble zmin, zmax, xc, yc, r, phi, mp, mn, dx, off, dh;
    guint *counts = NULL;
    guint i, mposp, mposn, res;
    const gdouble *d;
    gboolean valley;

    res = (guint)floor(0.8*sqrt(n) + 1.0);
    if (res < 8)
        return step1_init(xyz, n, param, estimcache);

    /* First we estimate the orientation (phi). */
    phi = estimate_projection_direction_red(xyz, n, estimcache);
    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    mean_line = gwy_data_line_new(res, 2.0*r, FALSE);
    counts = g_new(guint, res);
    gwy_data_line_set_offset(mean_line, -r);
    projection_to_line(xyz, n, phi, xc, yc, mean_line, NULL, counts);
    data_line_shorten(mean_line, counts, 4);
    d = gwy_data_line_get_data(mean_line);
    dx = gwy_data_line_get_dx(mean_line);
    off = gwy_data_line_get_offset(mean_line);

    mposp = mposn = res/2;
    mp = mn = 0.0;
    for (i = 1; i < res; i++) {
        if (!counts[i] || !counts[i-1])
           continue;
        dh = d[i] - d[i-1];
        if (dh > mn) {
            mn = dh;
            mposn = i;
        }
        else if (dh < mp) {
            mp = dh;
            mposp = i;
        }
    }

    valley = (mposp > mposn);
    param[0] = (valley ? -1 : 1)*(zmax - zmin);
    param[1] = dx*(MAX(mposp, mposn) - MIN(mposp, mposn) + 1);
    param[2] = 0.01*r;
    param[3] = valley ? zmax : zmin;
    param[4] = xc*cos(phi) - yc*sin(phi) + 0.5*(mposp + mposn)*dx + off;
    param[5] = gwy_canonicalize_angle(phi, FALSE, FALSE);
    param[6] = param[7] = 0.0;

    return TRUE;
}

/*******************************************************************************************************************
 *
 * Sphere
 *
 *******************************************************************************************************************/

static gdouble
sphere_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble kappa = param[3];
    gdouble r2k, t, val;

    x -= xc;
    y -= yc;
    /* Rewrite R - sqrt(R² - r²) as κ*r²/(1 + sqrt(1 - κ²r²)) where r² = x² + y² and κR = 1 to get nice behaviour in
     * the close-to-denegerate cases, including completely flat surface.  The expression 1.0/kappa is safe because we
     * cannot get to this branch for κ → 0 unless simultaneously r → ∞. */
    r2k = kappa*(x*x + y*y);
    t = 1.0 - kappa*r2k;
    if (t > 0.0)
        val = z0 + r2k/(1.0 + sqrt(t));
    else
        val = z0 + 2.0/kappa;

    return val;
}

static gboolean
sphere_init(const GwyXYZ *xyz, guint n, gdouble *param,
            FitShapeEstimateCache *estimcache)
{
    gdouble xc, yc, r, zmin, zmax, zmean;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);
    stat_z(xyz, n, &zmean, NULL, NULL, estimcache);

    param[0] = xc;
    param[1] = yc;
    if (fabs(zmean - zmin) > fabs(zmean - zmax)) {
        param[2] = zmax;
        param[3] = 2.0*(zmin - zmax)/(r*r);
    }
    else {
        param[2] = zmin;
        param[3] = 2.0*(zmax - zmin)/(r*r);
    }

    return TRUE;
}

/* Fit the data with a rotationally symmetric parabola and use its parameters for the spherical surface estimate. */
static gboolean
sphere_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                FitShapeEstimateCache *estimcache)
{
    gdouble xc, yc;
    /* Linear fit with functions 1, x, y and x²+y². */
    gdouble a[10], b[4];
    guint i;

    /* XXX: Handle the surrounding flat area, which can be a part of the function, better? */

    /* Using centered coodinates improves the condition number. */
    mean_x_y(xyz, n, &xc, &yc, estimcache);
    gwy_clear(a, 10);
    gwy_clear(b, 4);
    for (i = 0; i < n; i++) {
        gdouble x = xyz[i].x - xc, y = xyz[i].y - yc;
        gdouble r2 = x*x + y*y;

        b[0] += xyz[i].z;
        b[1] += x*xyz[i].z;
        b[2] += y*xyz[i].z;
        b[3] += r2*xyz[i].z;

        a[2] += x*x;
        a[4] += x*y;
        a[5] += y*y;
        a[6] += r2;
        a[7] += x*r2;
        a[8] += y*r2;
        a[9] += r2*r2;
    }
    a[0] = n;
    a[1] = a[3] = 0.0;

    param[0] = xc;
    param[1] = yc;
    param[2] = b[0]/n;
    param[3] = 0.0;

    if (!gwy_math_choleski_decompose(4, a))
        return FALSE;

    gwy_math_choleski_solve(4, a, b);

    param[3] = 2.0*b[3];
    if (param[3]) {
        param[0] = xc - b[1]/param[3];
        param[1] = yc - b[2]/param[3];
        param[2] = b[0] - 0.5*(b[1]*b[1] + b[2]*b[2])/param[3];
    }

    return TRUE;
}

static gdouble
sphere_calc_R(const gdouble *param)
{
    return 1.0/param[3];
}

static gdouble
sphere_calc_err_R(const gdouble *param,
                  const gdouble *param_err,
                  G_GNUC_UNUSED const gdouble *correl)
{
    return param_err[3]/(param[3]*param[3]);
}

static gdouble
sphere_calc_zc(const gdouble *param)
{
    return param[2] + 1.0/param[3];
}

static gdouble
sphere_calc_err_zc(const gdouble *param,
                   const gdouble *param_err,
                   const gdouble *correl)
{
    gdouble diff[G_N_ELEMENTS(sphere_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[2] = 1.0;
    diff[3] = -1.0/(param[3]*param[3]);
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

/*******************************************************************************************************************
 *
 * Half-sphere
 *
 *******************************************************************************************************************/

static gdouble
hsphere_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble kappa = param[3];
    gdouble r2k, t, val;

    x -= xc;
    y -= yc;
    /* Rewrite R - sqrt(R² - r²) as κ*r²/(1 + sqrt(1 - κ²r²)) where r² = x² + y² and κR = 1 to get nice behaviour in
     * the close-to-denegerate cases, including completely flat surface.  The expression 1.0/kappa is safe because we
     * cannot get to this branch for κ → 0 unless simultaneously r → ∞. */
    r2k = kappa*(x*x + y*y);
    t = 1.0 - kappa*r2k;
    if (t > 0.0)
        val = z0 + r2k/(1.0 + sqrt(t));
    else
        val = z0 + 1.0/kappa;

    return val;
}

static gboolean
hsphere_init(const GwyXYZ *xyz, guint n, gdouble *param,
             FitShapeEstimateCache *estimcache)
{
    return sphere_init(xyz, n, param, estimcache);
}

/* Fit the data with a rotationally symmetric parabola and use its parameters for the spherical surface estimate. */
static gboolean
hsphere_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                 FitShapeEstimateCache *estimcache)
{
    return sphere_estimate(xyz, n, param, estimcache);
}

static gdouble
hsphere_calc_R(const gdouble *param)
{
    return 1.0/param[3];
}

static gdouble
hsphere_calc_err_R(const gdouble *param,
                   const gdouble *param_err,
                   G_GNUC_UNUSED const gdouble *correl)
{
    return param_err[3]/(param[3]*param[3]);
}

static gdouble
hsphere_calc_zc(const gdouble *param)
{
    return param[2] + 1.0/param[3];
}

static gdouble
hsphere_calc_err_zc(const gdouble *param,
                    const gdouble *param_err,
                    const gdouble *correl)
{
    gdouble diff[G_N_ELEMENTS(hsphere_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[2] = 1.0;
    diff[3] = -1.0/(param[3]*param[3]);
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

/*******************************************************************************************************************
 *
 * Cylinder (lying)
 *
 *******************************************************************************************************************/

/* XXX: We might want
 * - finite-size cylinders
 * - cylinders cut elsewhere than in the half
 * Unfortunately, the derivatives by the corresponding parameters are often zero (or something odd) because of the
 * sharp boundaries in the function that either catch a data point inside some region or not. */
static gdouble
hcylinder_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble x0 = param[0];
    gdouble z0 = param[1];
    gdouble kappa = param[2];
    gdouble phi = param[3];
    gdouble bparallel = param[4];
    gdouble r2k, t, val, cphi, sphi;

    _gwy_sincos(phi, &sphi, &cphi);
    z0 += bparallel*(x*sphi + y*cphi);
    t = x*cphi - y*sphi - x0;
    r2k = kappa*t*t;
    t = 1.0 - kappa*r2k;
    if (t > 0.0)
        val = z0 + r2k/(1.0 + sqrt(t));
    else
        val = z0 + 2.0/kappa;

    return val;
}

static gboolean
hcylinder_init(const GwyXYZ *xyz, guint n, gdouble *param,
               FitShapeEstimateCache *estimcache)
{
    gdouble xc, yc, r, zmin, zmax, zmean;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);
    stat_z(xyz, n, &zmean, NULL, NULL, estimcache);

    param[0] = xc;
    if (fabs(zmean - zmin) > fabs(zmean - zmax)) {
        param[1] = zmax;
        param[2] = 2.0*(zmin - zmax)/(r*r);
    }
    else {
        param[1] = zmin;
        param[2] = 2.0*(zmax - zmin)/(r*r);
    }
    param[3] = 0.0;
    param[4] = 0.0;

    return TRUE;
}

/* Fit the data with a rotationally symmetric parabola and use its parameters for the spherical surface estimate. */
static gboolean
hcylinder_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                   FitShapeEstimateCache *estimcache)
{
    GwyDataLine *mean_line;
    gdouble xc, yc, r, zmin, zmax, t, phi;
    guint *counts = NULL;
    guint i, mpos, res;
    gint d2sign;
    const gdouble *d;

    /* First we estimate the orientation (phi). */
    phi = estimate_projection_direction_red(xyz, n, estimcache);

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    res = (guint)floor(0.8*sqrt(n) + 1.0);
    mean_line = gwy_data_line_new(res, 2.0*r, FALSE);
    counts = g_new(guint, res);
    gwy_data_line_set_offset(mean_line, -r);
    projection_to_line(xyz, n, phi, xc, yc, mean_line, NULL, counts);
    d = gwy_data_line_get_data(mean_line);

    d2sign = 0;
    for (i = 1; i < res-1; i++) {
        if (d[i] < 0.5*(d[i-1] + d[i+1]))
            d2sign += 1;
        else if (d[i] > 0.5*(d[i-1] + d[i+1]))
            d2sign -= 1;
    }
    gwy_debug("d2sign %d", d2sign);

    if (d2sign > 0)
        gwy_data_line_multiply(mean_line, -1.0);
    mpos = 0;
    for (i = 1; i < res; i++) {
        if (d[i] > d[mpos])
            mpos = i;
    }

    t = 2.0*r/res*mpos - r;
    param[0] = ((xc + t*cos(phi))*cos(phi) - (yc - t*sin(phi))*sin(phi));

    if (d2sign > 0) {
        param[1] = zmin;
        param[2] = 4.0*(zmax - zmin)/(r*r);
    }
    else {
        param[1] = zmax;
        param[2] = 4.0*(zmin - zmax)/(r*r);
    }
    param[3] = phi;
    param[4] = 0.0;

    g_free(counts);
    g_object_unref(mean_line);

    return TRUE;
}

static gdouble
hcylinder_calc_R(const gdouble *param)
{
    return 1.0/param[2];
}

static gdouble
hcylinder_calc_err_R(const gdouble *param,
                     const gdouble *param_err,
                     G_GNUC_UNUSED const gdouble *correl)
{
    return param_err[2]/(param[2]*param[2]);
}

/*******************************************************************************************************************
 *
 * Grating (simple)
 *
 *******************************************************************************************************************/

static gdouble
grating_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble L = fabs(param[0]);
    gdouble h = param[1];
    gdouble p = fabs(param[2]);
    gdouble z0 = param[3];
    gdouble x0 = param[4];
    gdouble phi = param[5];
    gdouble c = param[6];
    gdouble t, Lp2, val, cphi, sphi;

    Lp2 = 0.5*L*p;
    if (G_UNLIKELY(!Lp2))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi - x0 + Lp2;
    t = (t - L*floor(t/L))/Lp2 - 1.0;
    if (fabs(t) < 1.0)
        val = z0 + h*(1.0 - gwy_coshm1(c*t)/gwy_coshm1(c));
    else
        val = z0;

    return val;
}

static gboolean
grating_init(const GwyXYZ *xyz, guint n, gdouble *param,
             FitShapeEstimateCache *estimcache)
{
    gdouble xc, yc, r, zmin, zmax;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    param[0] = r/4.0;
    param[1] = zmax - zmin;
    param[2] = 0.5;
    param[3] = zmin;
    param[4] = 0.0;
    param[5] = 0.0;
    param[6] = 5.0;

    return TRUE;
}

static gboolean
grating_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                 FitShapeEstimateCache *estimcache)
{
    gdouble t;

    /* Just initialise the percentage and shape with some sane defaults. */
    param[2] = 0.5;
    param[6] = 5.0;

    /* Simple height parameter estimate. */
    range_z(xyz, n, param+3, param+1, estimcache);
    t = param[1] - param[3];
    param[1] = 0.9*t;
    param[3] += 0.05*t;

    /* First we estimate the orientation (phi). */
    param[5] = estimate_projection_direction_red(xyz, n, estimcache);

    /* Then we extract a representative profile with this orientation. */
    return estimate_period_and_phase(xyz, n, param[5], param + 0, param + 4, estimcache);
}

/*******************************************************************************************************************
 *
 * Grating (3-level)
 *
 *******************************************************************************************************************/

static gdouble
grating3_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble L = fabs(param[0]);
    gdouble h1 = fabs(param[1]);
    gdouble h2 = fabs(param[2]);
    gdouble h3 = fabs(param[3]);
    gdouble p = fmin(fabs(param[4]), 1.0);
    gdouble q1 = fabs(param[5]);
    gdouble q2 = fabs(param[6]);
    gdouble q3 = fabs(param[7]);
    gdouble z0 = param[8];
    gdouble x0 = param[9];
    gdouble phi = param[10];
    gdouble t, Lp2, cphi, sphi, Ll, Lu;

    Lp2 = 0.5*L*p;
    if (G_UNLIKELY(!Lp2))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi - x0 + Lp2;
    t -= L*floor(t/L) + Lp2;
    t = fabs(t);

    Lu = Lp2;
    if (t >= Lu)
        return z0;

    Ll = Lu;
    Lu = Ll/(1.0 + q1);
    if (t >= Lu) {
        if (G_UNLIKELY(Lu == Ll))
            return z0;
        return z0 + h1/(Lu - Ll)*(t - Ll);
    }

    Ll = Lu;
    Lu = Ll/(1.0 + q2);
    z0 += h1;
    if (t >= Lu) {
        if (G_UNLIKELY(Lu == Ll))
            return z0;
        return z0 + h2/(Lu - Ll)*(t - Ll);
    }

    Ll = Lu;
    Lu = Ll/(1.0 + q3);
    z0 += h2;
    if (t >= Lu) {
        if (G_UNLIKELY(Lu == Ll))
            return z0;
        return z0 + h3/(Lu - Ll)*(t - Ll);
    }

    return z0 + h3;
}

static gboolean
grating3_init(const GwyXYZ *xyz, guint n, gdouble *param,
              FitShapeEstimateCache *estimcache)
{
    gdouble xc, yc, r, zmin, zmax;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    param[0] = r/4.0;
    param[1] = 0.1*(zmax - zmin);
    param[2] = 0.8*(zmax - zmin);
    param[3] = 0.1*(zmax - zmin);
    param[4] = 0.7;
    param[5] = 0.5;
    param[6] = 0.2;
    param[7] = 0.5;
    param[8] = zmin;
    param[9] = 0.0;
    param[10] = 0.0;

    return TRUE;
}

static gboolean
grating3_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                  FitShapeEstimateCache *estimcache)
{
    gdouble zmin, zmax;

    /* Just initialise the percentage and shape with some sane defaults. */
    param[4] = 0.7;
    param[5] = 0.5;
    param[6] = 0.2;
    param[7] = 0.5;

    /* Simple height parameter estimate. */
    range_z(xyz, n, &zmin, &zmax, estimcache);
    param[1] = 0.1*(zmax - zmin);
    param[2] = 0.8*(zmax - zmin);
    param[3] = 0.1*(zmax - zmin);
    param[8] = zmin;

    /* First we estimate the orientation (phi). */
    param[10] = estimate_projection_direction_red(xyz, n, estimcache);

    /* Then we extract a representative profile with this orientation. */
    return estimate_period_and_phase(xyz, n, param[10], param + 0, param + 9, estimcache);
}

static gdouble
grating3_calc_h(const gdouble *param)
{
    return param[1] + param[2] + param[3];
}

static gdouble
grating3_calc_err_h(G_GNUC_UNUSED const gdouble *param,
                    const gdouble *param_err,
                    const gdouble *correl)
{
    gdouble diff[G_N_ELEMENTS(grating3_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[1] = diff[2] = diff[3] = 1.0;
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
grating3_calc_L0(const gdouble *param)
{
    return param[0]*param[4];
}

static gdouble
grating3_calc_err_L0(const gdouble *param,
                     const gdouble *param_err,
                     const gdouble *correl)
{
    gdouble diff[G_N_ELEMENTS(grating3_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[0] = param[4];
    diff[4] = param[0];
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
grating3_calc_L1(const gdouble *param)
{
    return grating3_calc_L0(param)/(1.0 + param[5]);
}

static gdouble
grating3_calc_err_L1(const gdouble *param,
                     const gdouble *param_err,
                     const gdouble *correl)
{
    gdouble L1 = grating3_calc_L1(param);
    gdouble diff[G_N_ELEMENTS(grating3_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[0] = L1/param[0];
    diff[4] = L1/param[4];
    diff[5] = -L1/(1.0 + param[5]);
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
grating3_calc_L2(const gdouble *param)
{
    return grating3_calc_L1(param)/(1.0 + param[6]);
}

static gdouble
grating3_calc_err_L2(const gdouble *param,
                     const gdouble *param_err,
                     const gdouble *correl)
{
    gdouble L2 = grating3_calc_L2(param);
    gdouble diff[G_N_ELEMENTS(grating3_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[0] = L2/param[0];
    diff[4] = L2/param[4];
    diff[5] = -L2/(1.0 + param[5]);
    diff[6] = -L2/(1.0 + param[6]);
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
grating3_calc_L3(const gdouble *param)
{
    return grating3_calc_L2(param)/(1.0 + param[7]);
}

static gdouble
grating3_calc_err_L3(const gdouble *param,
                     const gdouble *param_err,
                     const gdouble *correl)
{
    gdouble L3 = grating3_calc_L3(param);
    gdouble diff[G_N_ELEMENTS(grating3_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[0] = L3/param[0];
    diff[4] = L3/param[4];
    diff[5] = -L3/(1.0 + param[5]);
    diff[6] = -L3/(1.0 + param[6]);
    diff[7] = -L3/(1.0 + param[7]);
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

/*******************************************************************************************************************
 *
 * Holes
 *
 *******************************************************************************************************************/

static inline gdouble
hole_radial_intersection(gdouble q, gdouble A, gdouble R)
{
    gdouble A1q = A*(1.0 - q);
    gdouble q21 = 1.0 + q*q;
    gdouble D = R*R*q21 - A1q*A1q;
    gdouble sqrtD = sqrt(MAX(D, 0.0));
    gdouble x = (1.0 + q)*A + sqrtD;
    return x/sqrt(q21);
}

static gdouble
hole_shape(gdouble x, gdouble y, gdouble size, gdouble slope, gdouble roundness)
{
    gdouble rx, ry, r, rr, rsz;
    gdouble v = 0.0;

    if (roundness) {
        rsz = roundness*size;
        rx = fabs(x) - (size - rsz);
        ry = fabs(y) - (size - rsz);
        r = MAX(rx, ry);
        rr = MIN(rx, ry);
        if (r <= 0.0 || (r <= rsz && rr <= 0.0) || rx*rx + ry*ry <= rsz*rsz)
            v = -1.0;
        else if (slope) {
            gdouble ss = size + slope;
            rsz = roundness*ss;
            rx = fabs(x) - (ss - rsz);
            ry = fabs(y) - (ss - rsz);
            r = MAX(rx, ry);
            rr = MIN(rx, ry);
            if (r <= 0.0 || (r <= rsz && rr <= 0.0) || rx*rx + ry*ry <= rsz*rsz) {
                gdouble q = (rr + ss - rsz)/(r + ss - rsz);
                if (q <= 1.0 - roundness)
                    v = (r - rsz)/slope;
                else {
                    r = hole_radial_intersection(q, ss - rsz, rsz);
                    rr = hole_radial_intersection(q, size - roundness*size, roundness*size);
                    v = (sqrt(x*x + y*y) - r)/(r - rr);
                }
            }
        }
    }
    else {
        rx = fabs(x) - size;
        ry = fabs(y) - size;
        r = MAX(rx, ry);
        if (r <= 0.0)
            v = -1.0;
        else if (r < slope)
            v = (r - slope)/slope;
    }

    return v;
}

static gdouble
holes_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble L = fabs(param[3]);
    gdouble p = param[4];
    gdouble h = param[5];
    gdouble s = param[6];
    gdouble r = param[7];
    gdouble phi = param[8];
    gdouble t, cphi, sphi;

    x -= xc;
    y -= yc;

    if (G_UNLIKELY(L*p == 0.0))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    x -= L*floor(x/L) + 0.5*L;
    y -= L*floor(y/L) + 0.5*L;

    p = 1.0/(1.0 + fabs(p));
    /* Map zero values to no roundness and no slope. */
    s = fabs(s)/(1.0 + fabs(s));
    r = fabs(r)/(1.0 + fabs(r));

    return z0 + h*hole_shape(fabs(x), fabs(y), (1.0 - s)*0.5*L*p, s*0.5*L*p, r);
}

static gboolean
holes_init(const GwyXYZ *xyz, guint n, gdouble *param,
           FitShapeEstimateCache *estimcache)
{
    gdouble xc, yc, r, zmin, zmax;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    param[0] = 0;
    param[1] = 0;
    param[2] = zmin;
    param[3] = r/4.0;
    param[4] = 1.0;
    param[5] = zmax - zmin;
    param[6] = 0.1;
    param[7] = 0.1;
    param[8] = 0.0;

    return TRUE;
}

static gboolean
holes_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
               FitShapeEstimateCache *estimcache)
{
    GwyDataLine *mean_line, *rms_line;
    gdouble xc, yc, r, phi, L1, L2, zmin, zmax, u, v;
    guint *counts;
    guint res;
    gboolean ok1, ok2;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    param[4] = 1.0;
    param[6] = 0.1;
    param[7] = 0.1;

    param[8] = phi = estimate_projection_direction_red(xyz, n, estimcache);

    ok1 = estimate_period_and_phase(xyz, n, phi, &L1, &u, estimcache);
    ok2 = estimate_period_and_phase(xyz, n, phi - 0.5*G_PI, &L2, &v, estimcache);
    param[3] = 0.5*(L1 + L2);
    param[0] = u*cos(phi) + v*sin(phi);
    param[1] = u*sin(phi) + v*cos(phi);

    /* Estimate h sign: do projection and if smaller values correlate with
     * large rms then it is holes (h > 0); if larger values correlate with
     * large rms then it is bumps (h < 0).  */
    res = (guint)floor(0.8*sqrt(n) + 1.0);
    mean_line = gwy_data_line_new(res, 2.0*r, FALSE);
    gwy_data_line_set_offset(mean_line, -r);
    rms_line = gwy_data_line_new_alike(mean_line, FALSE);
    counts = g_new(guint, res);
    projection_to_line(xyz, n, phi, xc, yc, mean_line, rms_line, counts);
    data_line_shorten(mean_line, counts, 4);
    data_line_shorten(rms_line, counts, 4);
    g_free(counts);

    if (data_line_pearson_coeff(mean_line, rms_line) <= 0.0) {
        gwy_debug("holes");
        param[2] = zmax;
        param[5] = zmax - zmin;
    }
    else {
        gwy_debug("bumps");
        param[2] = zmin;
        param[5] = zmin - zmax;
        /* estimate_period_and_phase() finds maxima but we want minima here */
        param[0] += 0.5*param[3];
        param[1] += 0.5*param[3];
    }

    g_object_unref(rms_line);
    g_object_unref(mean_line);

    return ok1 && ok2;
}

static gdouble
holes_calc_wouter(const gdouble *param)
{
    return param[3]/(1.0 + fabs(param[4]));
}

static gdouble
holes_calc_err_wouter(const gdouble *param,
                      const gdouble *param_err,
                      const gdouble *correl)
{
    gdouble diff[G_N_ELEMENTS(holes_params)];
    gdouble wouter = holes_calc_wouter(param);

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[3] = wouter/param[3];
    diff[4] = -wouter/(1.0 + fabs(param[4]));
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
holes_calc_winner(const gdouble *param)
{
    return param[3]/(1.0 + fabs(param[4]))/(1.0 + fabs(param[6]));
}

static gdouble
holes_calc_err_winner(const gdouble *param,
                      const gdouble *param_err,
                      const gdouble *correl)
{
    gdouble diff[G_N_ELEMENTS(holes_params)];
    gdouble winner = holes_calc_winner(param);

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[3] = winner/param[3];
    diff[4] = -winner/(1.0 + fabs(param[4]));
    diff[6] = -winner/(1.0 + fabs(param[6]));
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
holes_calc_R(const gdouble *param)
{
    return 0.5*param[3]/(1.0 + fabs(param[4]))*fabs(param[7])/(1.0 + fabs(param[7]));
}

static gdouble
holes_calc_err_R(const gdouble *param,
                 const gdouble *param_err,
                 const gdouble *correl)
{
    gdouble diff[G_N_ELEMENTS(holes_params)];
    gdouble R = holes_calc_R(param);
    gdouble t = 1.0 + fabs(param[7]), u = 1.0 + fabs(param[4]);

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[3] = R/param[3];
    diff[4] = -R/u;
    /* Do this directly to avoid division by possibly zero r. */
    diff[7] = 0.5*param[3]/(u*t*t);
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

/*******************************************************************************************************************
 *
 * Ring
 *
 *******************************************************************************************************************/

static gdouble
pring_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble R = param[3];
    gdouble w = fabs(param[4]);
    gdouble h = param[5];
    gdouble s = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble r, r2, s_h, rinner, router;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;
    r2 = x*x + y*y;

    if (G_UNLIKELY(w == 0.0))
        return r2 <= R*R ? z0 - 0.5*s : z0 + 0.5*s;

    if (G_UNLIKELY(h == 0.0)) {
        if (r2 >= R*R)
            return z0 + 0.5*s;
        else if (r2 < (R - w)*(R - w))
            return z0 - 0.5*s;

        r = (R - sqrt(r2))/w;
        return z0 + s*(0.5 - r*r);
    }

    r = sqrt(r2) - R;
    s_h = s/h;
    rinner = sqrt(fmax(1.0 + 0.5*s_h, 0.0));
    router = sqrt(fmax(1.0 - 0.5*s_h, 0.0));

    rinner *= -0.5*w;
    if (r <= rinner)
        return z0 - 0.5*s;

    router *= 0.5*w;
    if (r >= router)
        return z0 + 0.5*s;

    r *= 2.0/w;
    return z0 + h*(1.0 - r*r);
}

static gboolean
pring_init(const GwyXYZ *xyz, guint n, gdouble *param,
           FitShapeEstimateCache *estimcache)
{
    gdouble xc, yc, r, zmin, zmax;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    param[0] = xc;
    param[1] = yc;
    param[2] = zmin;
    param[3] = r/3.0;
    param[4] = r/12.0;
    param[5] = zmax - zmin;
    param[6] = (zmax - zmin)/12.0;
    param[7] = 0.0;
    param[8] = 0.0;

    return TRUE;
}

static gboolean
pring_estimate_projection(const GwyXYZ *xyz, guint n,
                          gdouble xc, gdouble yc, gdouble r,
                          gboolean vertical, gboolean upwards,
                          GwyDataLine *proj, GwyXY *projdata, guint *counts,
                          GwyPeaks *peaks, gdouble *param)
{
    guint i, res;
    const gdouble *d;
    gdouble c, real, off, width[2];

    c = (vertical ? yc : xc);
    gwy_data_line_set_real(proj, 2.0*r);
    gwy_data_line_set_offset(proj, -r);
    projection_to_line(xyz, n, vertical ? -0.5*G_PI : 0.0, xc, yc, proj, NULL, counts);
    data_line_shorten(proj, counts, 4);

    if (!upwards)
        gwy_data_line_multiply(proj, -1.0);

    res = gwy_data_line_get_res(proj);
    real = gwy_data_line_get_real(proj);
    off = gwy_data_line_get_offset(proj);
    d = gwy_data_line_get_data(proj);
    for (i = 0; i < res; i++) {
        projdata[i].x = c + off + (i + 0.5)*real/res;
        projdata[i].y = d[i];
    }
    if (gwy_peaks_analyze_xy(peaks, projdata, res, 2) != 2)
        return FALSE;

    gwy_peaks_get_quantity(peaks, GWY_PEAK_ABSCISSA, param);
    gwy_peaks_get_quantity(peaks, GWY_PEAK_WIDTH, width);
    param[2] = 0.5*(width[0] + width[1]);
    return TRUE;
}

static gboolean
pring_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
               FitShapeEstimateCache *estimcache)
{
    GwyDataLine *proj;
    GwyXY *projdata;
    GwyPeaks *peaks;
    gdouble xc, yc, r, zmin, zmax, zskew;
    gdouble xestim[3], yestim[3];
    guint *counts;
    gboolean ok1, ok2;
    guint res;

    circumscribe_x_y(xyz, n, &xc, &yc, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);
    stat_z(xyz, n, NULL, NULL, &zskew, estimcache);
    res = (guint)floor(0.8*sqrt(n) + 1.0);
    if (zskew < 0.0)
        GWY_SWAP(gdouble, zmin, zmax);

    proj = gwy_data_line_new(res, 2.0*r, FALSE);
    counts = g_new(guint, res);
    projdata = g_new(GwyXY, res);

    peaks = gwy_peaks_new();
    gwy_peaks_set_order(peaks, GWY_PEAK_ORDER_ABSCISSA);
    gwy_peaks_set_background(peaks, GWY_PEAK_BACKGROUND_MMSTEP);

    gwy_data_line_resample(proj, res, GWY_INTERPOLATION_NONE);
    ok1 = pring_estimate_projection(xyz, n, xc, yc, r, FALSE, zskew >= 0.0, proj, projdata, counts, peaks, xestim);

    gwy_data_line_resample(proj, res, GWY_INTERPOLATION_NONE);
    ok2 = pring_estimate_projection(xyz, n, xc, yc, r, TRUE, zskew >= 0.0, proj, projdata, counts, peaks, yestim);

    g_free(counts);
    g_object_unref(proj);
    gwy_peaks_free(peaks);

    if (!ok1 || !ok2)
        return FALSE;

    param[0] = 0.5*(xestim[0] + xestim[1]);
    param[1] = 0.5*(yestim[0] + yestim[1]);
    param[2] = zmin;
    param[3] = 0.25*(yestim[1] - yestim[0] + xestim[1] - xestim[0]);
    /* A bit too high value is OK because at least the estimated function does not miss the ring completely. */
    param[4] = 1.5*(xestim[2] + yestim[2]);
    param[5] = zmax - zmin;
    param[6] = (zmax - zmin)/12.0;
    param[7] = 0.0;
    param[8] = 0.0;

    return TRUE;
}

/*******************************************************************************************************************
 *
 * Gaussian
 *
 *******************************************************************************************************************/

static gdouble
gaussian_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble sigma = param[4];
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble t, val, cphi, sphi, s2;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    s2 = sigma*sigma;
    if (G_UNLIKELY(!s2 || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    t = 0.5*(x*x*a + y*y/a)/s2;
    val = z0 + h*exp(-t);

    return val;
}

static void
gaussian_diff(gdouble x, gdouble y,
              const gdouble *param,
              G_GNUC_UNUSED const gboolean *fixed_param,
              gdouble *der)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    /* gdouble z0 = param[2]; */
    gdouble h = param[3];
    gdouble sigma = fabs(param[4]);
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble xx, yy, cphi, sphi, aa, s2, t2, g;

    s2 = sigma*sigma;
    if (G_UNLIKELY(!s2 || !a)) {
        gwy_clear(der, G_N_ELEMENTS(gaussian_params));
        return;
    }

    x -= xc;
    y -= yc;

    der[2] = 1;
    der[7] = x;
    der[8] = y;

    /* We need to calculate these for any of the following parameters.  And once we have them, the rest is just simple
     * algebra.  So do not bother checking whether parameters are fixed or not. */
    _gwy_sincos(phi, &sphi, &cphi);
    xx = x*cphi - y*sphi;
    yy = x*sphi + y*cphi;
    t2 = xx*xx*a + yy*yy/a;
    aa = a - 1.0/a;
    g = exp(-0.5*t2/s2);

    der[0] = -bx + h*g/s2 * (a*xx*cphi + yy/a*sphi);
    der[1] = -by - h*g/s2 * (a*xx*sphi - yy/a*cphi);
    der[3] = g;
    der[4] = h*g * t2/(s2*sigma);
    der[5] = 0.5*h*g/s2 * (yy*yy/(a*a) - xx*xx);
    der[6] = h*g/s2 * aa*xx*yy;
}

static gboolean
gaussian_init(const GwyXYZ *xyz, guint n, gdouble *param,
              FitShapeEstimateCache *estimcache)
{
    param[7] = 0.0;
    param[8] = 0.0;
    return common_bump_feature_init(xyz, n,
                                    param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                    estimcache);
}

static gboolean
gaussian_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                  FitShapeEstimateCache *estimcache)
{
    param[7] = 0.0;
    param[8] = 0.0;
    return common_bump_feature_estimate(xyz, n,
                                        param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                        estimcache);
}

static gdouble
gaussian_calc_sigma1(const gdouble *param)
{
    return param[4]/sqrt(fabs(param[5]));
}

static gdouble
gaussian_calc_err_sigma1(const gdouble *param,
                         const gdouble *param_err,
                         const gdouble *correl)
{
    gdouble sigma1 = gaussian_calc_sigma1(param);
    gdouble diff[G_N_ELEMENTS(gaussian_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[4] = sigma1/param[4];
    diff[5] = -0.5*sigma1/param[5];
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
gaussian_calc_sigma2(const gdouble *param)
{
    return param[4]*sqrt(fabs(param[5]));
}

static gdouble
gaussian_calc_err_sigma2(const gdouble *param,
                         const gdouble *param_err,
                         const gdouble *correl)
{
    gdouble sigma1 = gaussian_calc_sigma1(param);
    gdouble diff[G_N_ELEMENTS(gaussian_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[4] = sigma1/param[4];
    diff[5] = 0.5*sigma1/param[5];
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

/*******************************************************************************************************************
 *
 * Lorentzian
 *
 *******************************************************************************************************************/

static gdouble
lorentzian_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble b = param[4];
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble t, val, cphi, sphi, b2;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    b2 = b*b;
    if (G_UNLIKELY(!b2 || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    t = (x*x*a + y*y/a)/b2;
    val = z0 + h/(1.0 + t);

    return val;
}

static void
lorentzian_diff(gdouble x, gdouble y,
                const gdouble *param,
                G_GNUC_UNUSED const gboolean *fixed_param,
                gdouble *der)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    /* gdouble z0 = param[2]; */
    gdouble h = param[3];
    gdouble b = fabs(param[4]);
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble xx, yy, cphi, sphi, aa, b2, t2, L;

    b2 = b*b;
    if (G_UNLIKELY(!b2 || !a)) {
        gwy_clear(der, G_N_ELEMENTS(lorentzian_params));
        return;
    }

    x -= xc;
    y -= yc;

    der[2] = 1;
    der[7] = x;
    der[8] = y;

    /* We need to calculate these for any of the following parameters.  And once we have them, the rest is just simple
     * algebra.  So do not bother checking whether parameters are fixed or not. */
    _gwy_sincos(phi, &sphi, &cphi);
    xx = x*cphi - y*sphi;
    yy = x*sphi + y*cphi;
    t2 = (xx*xx*a + yy*yy/a)/b2;
    aa = a - 1.0/a;
    L = 1.0/(1.0 + t2);

    der[0] = -bx + 2.0*h*L*L/b2 * (a*xx*cphi + yy/a*sphi);
    der[1] = -by - 2.0*h*L*L/b2 * (a*xx*sphi - yy/a*cphi);
    der[3] = L;
    der[4] = h*L*L/b2 * 2.0*b*t2;
    der[5] = h*L*L/b2 * (yy*yy/(a*a) - xx*xx);
    der[6] = h*L*L/b2 * 2.0*aa*xx*yy;
}

static gboolean
lorentzian_init(const GwyXYZ *xyz, guint n, gdouble *param,
                FitShapeEstimateCache *estimcache)
{
    param[7] = 0.0;
    param[8] = 0.0;
    return common_bump_feature_init(xyz, n,
                                    param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                    estimcache);
}

static gboolean
lorentzian_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                    FitShapeEstimateCache *estimcache)
{
    param[7] = 0.0;
    param[8] = 0.0;
    return common_bump_feature_estimate(xyz, n,
                                        param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                        estimcache);
}

static gdouble
lorentzian_calc_b1(const gdouble *param)
{
    return param[4]/sqrt(fabs(param[5]));
}

static gdouble
lorentzian_calc_err_b1(const gdouble *param,
                       const gdouble *param_err,
                       const gdouble *correl)
{
    gdouble b1 = lorentzian_calc_b1(param);
    gdouble diff[G_N_ELEMENTS(lorentzian_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[4] = b1/param[4];
    diff[5] = -0.5*b1/param[5];
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
lorentzian_calc_b2(const gdouble *param)
{
    return param[4]*sqrt(fabs(param[5]));
}

static gdouble
lorentzian_calc_err_b2(const gdouble *param,
                       const gdouble *param_err,
                       const gdouble *correl)
{
    gdouble b1 = lorentzian_calc_b1(param);
    gdouble diff[G_N_ELEMENTS(lorentzian_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[4] = b1/param[4];
    diff[5] = 0.5*b1/param[5];
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

/*******************************************************************************************************************
 *
 * Exponential
 *
 *******************************************************************************************************************/

static gdouble
exponential_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble b = param[4];
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble t, val, cphi, sphi, b2;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    b2 = b*b;
    if (G_UNLIKELY(!b2 || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    t = (x*x*a + y*y/a)/b2;
    val = z0 + h*exp(-sqrt(t));

    return val;
}

static void
exponential_diff(gdouble x, gdouble y,
                 const gdouble *param,
                 G_GNUC_UNUSED const gboolean *fixed_param,
                 gdouble *der)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    /* gdouble z0 = param[2]; */
    gdouble h = param[3];
    gdouble b = fabs(param[4]);
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble xx, yy, cphi, sphi, aa, t, E;

    if (G_UNLIKELY(!b || !a)) {
        gwy_clear(der, G_N_ELEMENTS(exponential_params));
        return;
    }

    x -= xc;
    y -= yc;

    der[2] = 1;
    der[7] = x;
    der[8] = y;

    /* We need to calculate these for any of the following parameters.  And once we have them, the rest is just simple
     * algebra.  So do not bother checking whether parameters are fixed or not. */
    _gwy_sincos(phi, &sphi, &cphi);
    xx = x*cphi - y*sphi;
    yy = x*sphi + y*cphi;
    t = sqrt(xx*xx*a + yy*yy/a);
    /* The function is not differentiable at t = 0.  Just skip the point when calculating derivatives. */
    if (G_UNLIKELY(!t)) {
        gwy_clear(der, G_N_ELEMENTS(exponential_params));
        return;
    }
    aa = a - 1.0/a;
    E = exp(-t/b);

    der[0] = -bx + h*E/(b*t) * (a*xx*cphi + yy/a*sphi);
    der[1] = -by - h*E/(b*t) * (a*xx*sphi - yy/a*cphi);
    der[3] = E;
    der[4] = h*E*t/(b*b);
    der[5] = h*E/(b*t) * 0.5*(yy*yy/(a*a) - xx*xx);
    der[6] = h*E/(b*t) * aa*xx*yy;
}

static gboolean
exponential_init(const GwyXYZ *xyz, guint n, gdouble *param,
                 FitShapeEstimateCache *estimcache)
{
    param[7] = 0.0;
    param[8] = 0.0;
    return common_bump_feature_init(xyz, n,
                                    param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                    estimcache);
}

static gboolean
exponential_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                     FitShapeEstimateCache *estimcache)
{
    param[7] = 0.0;
    param[8] = 0.0;
    return common_bump_feature_estimate(xyz, n,
                                        param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                        estimcache);
}

static gdouble
exponential_calc_b1(const gdouble *param)
{
    return param[4]/sqrt(fabs(param[5]));
}

static gdouble
exponential_calc_err_b1(const gdouble *param,
                        const gdouble *param_err,
                        const gdouble *correl)
{
    gdouble b1 = exponential_calc_b1(param);
    gdouble diff[G_N_ELEMENTS(exponential_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[4] = b1/param[4];
    diff[5] = -0.5*b1/param[5];
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
exponential_calc_b2(const gdouble *param)
{
    return param[4]*sqrt(fabs(param[5]));
}

static gdouble
exponential_calc_err_b2(const gdouble *param,
                        const gdouble *param_err,
                        const gdouble *correl)
{
    gdouble b1 = exponential_calc_b1(param);
    gdouble diff[G_N_ELEMENTS(exponential_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[4] = b1/param[4];
    diff[5] = 0.5*b1/param[5];
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

/*******************************************************************************************************************
 *
 * Pyramid (diamond)
 *
 *******************************************************************************************************************/

static gdouble
pyramidx_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble L = param[4];
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble t, val, cphi, sphi, q;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    if (G_UNLIKELY(!L || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    q = 0.5*L*sqrt(1.0 + a*a);
    x /= q;
    y *= a/q;
    t = fabs(x) + fabs(y);
    if (t < 1.0)
        val = z0 + h*(1.0 - t);
    else
        val = z0;

    return val;
}

static gboolean
pyramidx_init(const GwyXYZ *xyz, guint n, gdouble *param,
              FitShapeEstimateCache *estimcache)
{
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_init(xyz, n,
                                  param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                  estimcache);
    param[4] *= 2.0;

    return ok;
}

static gboolean
pyramidx_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                  FitShapeEstimateCache *estimcache)
{
    gboolean ok;

    /* XXX: The pyramid has minimum projection when oriented along x and y axes.  But not very deep.  Can we use it to
     * estimate phi? */
    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_estimate(xyz, n,
                                      param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                      estimcache);
    param[4] *= 2.0;

    return ok;
}

/*******************************************************************************************************************
 *
 * Pyramid (rectangle)
 *
 *******************************************************************************************************************/

static gdouble
pyramid_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble L = param[4];
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble t, val, cphi, sphi, q;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    if (G_UNLIKELY(!L || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    q = 0.5/G_SQRT2*L*sqrt(1.0 + a*a);
    x /= q;
    y *= a/q;
    t = fmax(fabs(x), fabs(y));
    if (t < 1.0)
        val = z0 + h*(1.0 - t);
    else
        val = z0;

    return val;
}

static gboolean
pyramid_init(const GwyXYZ *xyz, guint n, gdouble *param,
             FitShapeEstimateCache *estimcache)
{
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_init(xyz, n,
                                  param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                  estimcache);
    param[4] *= G_SQRT2;

    return ok;
}

static gboolean
pyramid_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                 FitShapeEstimateCache *estimcache)
{
    gboolean ok;

    /* XXX: The pyramid has minimum projection when oriented along x and y axes.  But not very deep.  Can we use it to
     * estimate phi? */
    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_estimate(xyz, n,
                                      param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                      estimcache);
    param[4] *= G_SQRT2;

    return ok;
}

/*******************************************************************************************************************
 *
 * Pyramid (3-sided)
 *
 *******************************************************************************************************************/

static gdouble
pyramid3_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble L = param[4];
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble t1, t2, t3, t, val, cphi, sphi, q;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    if (G_UNLIKELY(!L || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    q = L/G_SQRT2*sqrt(1.0 + a*a);
    x /= q;
    y *= a/q;
    t1 = 2.0*GWY_SQRT3*x;
    t2 = -(GWY_SQRT3*x - 3.0*y);
    t3 = -(GWY_SQRT3*x + 3.0*y);
    t = fmax(fmax(t1, t2), t3);
    if (t < 1.0)
        val = z0 + h*(1.0 - t);
    else
        val = z0;

    return val;
}

static gboolean
pyramid3_init(const GwyXYZ *xyz, guint n, gdouble *param,
              FitShapeEstimateCache *estimcache)
{
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_init(xyz, n,
                                  param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                  estimcache);
    param[4] *= 2.0;

    return ok;
}

static gboolean
pyramid3_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                  FitShapeEstimateCache *estimcache)
{
    gboolean ok;

    /* XXX: The 3-pyramid projection has 6-fold symmetry and essentially no change in the projection while rotating.
     * */
    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_estimate(xyz, n,
                                      param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                      estimcache);
    param[4] *= 2.0;

    return ok;
}

/*******************************************************************************************************************
 *
 * Cone
 *
 *******************************************************************************************************************/

static gdouble
cone_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble R = param[4];
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble t, val, cphi, sphi, q;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    if (G_UNLIKELY(!R || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    q = R/G_SQRT2*sqrt(1.0 + a*a);
    x /= q;
    y *= a/q;
    t = sqrt(x*x + y*y);
    if (t < 1.0)
        val = z0 + h*(1.0 - t);
    else
        val = z0;

    return val;
}

static gboolean
cone_init(const GwyXYZ *xyz, guint n, gdouble *param,
          FitShapeEstimateCache *estimcache)
{
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_init(xyz, n,
                                  param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                  estimcache);
    param[4] *= G_SQRT2;

    return ok;
}

static gboolean
cone_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
              FitShapeEstimateCache *estimcache)
{
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_estimate(xyz, n,
                                      param + 0, param + 1, param + 2, param + 3, param + 4, param + 5, param + 6,
                                      estimcache);
    param[4] *= G_SQRT2;

    return ok;
}

/*******************************************************************************************************************
 *
 * Parabolic bump
 *
 *******************************************************************************************************************/

static gdouble
parbump_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble C = fabs(param[4]);
    gdouble a = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble t, val, cphi, sphi;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    if (G_UNLIKELY(!C || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    t = fabs(h) - 0.5*C*(x*x*a + y*y/a);
    if (h >= 0.0)
        val = (t > 0.0 ? z0 + t : z0);
    else
        val = (t > 0.0 ? z0 - t : z0);

    return val;
}

static gboolean
parbump_init(const GwyXYZ *xyz, guint n, gdouble *param,
             FitShapeEstimateCache *estimcache)
{
    gdouble w;
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_init(xyz, n,
                                  param + 0, param + 1, param + 2, param + 3, &w, param + 5, param + 6,
                                  estimcache);
    param[4] = 2.0*param[3]/(w*w);

    return ok && w >= 0.0;
}

static gboolean
parbump_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                 FitShapeEstimateCache *estimcache)
{
    gdouble w;
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_estimate(xyz, n,
                                      param + 0, param + 1, param + 2, param + 3, &w, param + 5, param + 6,
                                      estimcache);
    param[4] = 2.0*param[3]/(w*w);

    return ok && w >= 0.0;
}

/*******************************************************************************************************************
 *
 * Gaussian broadened cone
 *
 *******************************************************************************************************************/

/* Two-point Padé approximation for the universal function describing shape of Gaussian broadened cone.  Parameter
 * t is squared unitless distance from the centre. */
static inline double
gauss_smoothed_cone(gdouble t)
{
    gdouble s = sqrt(1.0 + t);
    gdouble n = 1.2533141373155 + t*(1.51210888400198 + t*(1.039435684177492 + t*(0.3985211963686722
                                                                                  + t*0.2494937570579042)));
    gdouble d = 1 + t*(1.206488332798031 + t*(0.7668496843527843 + t*(0.3259017222687021 + t*0.1764187274793531)));

    return s*n/d - 1.2533141373155;
}

static gdouble
smthcone_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble a = fabs(param[4]);
    gdouble R = fabs(param[5]);
    gdouble bx = param[6];
    gdouble by = param[7];
    gdouble p, rho2, val;

    x -= xc;
    y -= yc;
    rho2 = x*x + y*y;
    z0 += bx*x + by*y;

    if (G_UNLIKELY(!h || !a))
        return z0;

    /* XXX: If someone uses this with the base visible then the transition to the base should be probably also smooth.
     * This shape is actually just a cone with smooth top.  */

    if (G_UNLIKELY(!R))
        val = a*sqrt(rho2);
    else {
        p = R*a*GWY_SQRT_PI/G_SQRT2;
        val = gauss_smoothed_cone(2*rho2/(p*p));
        val *= 0.5*p*a;
    }

    if (val >= fabs(h))
        return z0;

    return h > 0.0 ? z0 + (h - val) : z0 - (h - val);
}

static gboolean
smthcone_init(const GwyXYZ *xyz, guint n, gdouble *param,
              FitShapeEstimateCache *estimcache)
{
    gdouble w, a_unused, phi_unused;
    gboolean ok;

    param[6] = 0.0;
    param[7] = 0.0;
    ok = common_bump_feature_init(xyz, n,
                                  param + 0, param + 1, param + 2, param + 3, &w, &a_unused, &phi_unused,
                                  estimcache);
    param[4] = 0.5*param[3]/w;
    param[5] = 0.2*w;

    return ok && w >= 0.0;
}

static gboolean
smthcone_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                  FitShapeEstimateCache *estimcache)
{
    gdouble w, a_unused, phi_unused;
    gboolean ok;

    param[6] = 0.0;
    param[7] = 0.0;
    ok = common_bump_feature_estimate(xyz, n,
                                      param + 0, param + 1, param + 2, param + 3, &w, &a_unused, &phi_unused,
                                      estimcache);
    param[4] = 0.5*param[3]/w;
    param[5] = 0.2*w;

    return ok && w >= 0.0;
}

/*******************************************************************************************************************
 *
 * Gaussian broadened 4-sided pyramid
 *
 *******************************************************************************************************************/

/* Exact formula along the u axis for v = 1 (i.e. at infinity).
 * We express it in terms of x, not u.  */
static inline gdouble
edgeatinf_pyramid(gdouble x)
{
    return (x*erf(x) - x) + exp(-x*x)/GWY_SQRT_PI;
}

/* Rational approximation along the axis v for u = 0. */
static inline gdouble
ratapprox_p4_u0(gdouble v)
{
    gdouble a = 2;
    gdouble b = -2.04551283690553;
    gdouble c = 1.31381602667408;
    gdouble C = -0.051558889057279;
    gdouble D = 0.0893893685714812;
    gdouble B = 1 + (b + c) - (C + D);
    gdouble v2 = v*v;

    return (a + b*v + c*v2)/(1 + B*v + C*v2 + D*v2*v2*v);
}

static inline gdouble
infapprox_p4(gdouble x, gdouble y)
{
    gdouble a = -0.376061521499245;
    gdouble b = -1.61222579045509;
    gdouble c = 0.326787561382881;
    gdouble d = 2.01942772652492;
    gdouble e = -0.4446887406374;
    gdouble f = 0.120470096230961;
    gdouble g = -1.19014793215417;
    gdouble h = 0.64604889949942;
    gdouble i = 0.569963272807616;
    gdouble x2 = x*x;
    gdouble y2 = y*y;
    gdouble r2 = x2 + y2;
    gdouble xx = x - y;
    gdouble u = sqrt(1 - exp(-pow(xx, 1.0/0.6)));
    gdouble v = sqrt(1 - exp(-y*y));
    gdouble fx = edgeatinf_pyramid(xx/G_SQRT2);
    gdouble fy = ratapprox_p4_u0(v);
    gdouble u2 = u*u;
    gdouble v2 = v*v;
    gdouble cbase = a * u*v*(1-v)*(1-v)*(1-u);
    gdouble denom = (1 + b*u + c*v + d*u2 + e*v2 + f*u*v + g*u2*u + h*v2*v + i*u2*v2);
    gdouble w = exp(-4*r2);
    gdouble c4phi = (x2*x2 - 6*x2*y2 + y2*y2)/(r2*r2 + 1e-16);
    gdouble s = (r2 + (c4phi - 3)/24*r2*r2 + (5 - 3*c4phi)/240*r2*r2*r2)/GWY_SQRT_PI;
    gdouble rv = fx*fy + cbase/denom + M_SQRT2*x - 2.0/GWY_SQRT_PI;

    return w*s + (1 - w)*rv;
}

static gdouble
smthpyr4_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble a = fabs(param[4]);
    gdouble R = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble sigmaq2, val, t, cphi, sphi;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    if (G_UNLIKELY(!h || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    /* fold coordinates to the unique octant of the plane */
    y = fabs(y);
    x = fabs(x);
    GWY_ORDER(gdouble, x, y);

    if (G_UNLIKELY(!R)) {
        val = x*a;
    }
    else {
        sigmaq2 = G_SQRT2/GWY_SQRT_PI*R*a;
        val = infapprox_p4(x/sigmaq2, y/sigmaq2);
        val *= sigmaq2*a/G_SQRT2;
    }

    if (val >= fabs(h))
        return z0;

    return h > 0.0 ? z0 + (h - val) : z0 - (h - val);
}

static gboolean
smthpyr4_init(const GwyXYZ *xyz, guint n, gdouble *param,
              FitShapeEstimateCache *estimcache)
{
    gdouble w, a_unused;
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_init(xyz, n,
                                  param + 0, param + 1, param + 2, param + 3, &w, &a_unused, param + 6,
                                  estimcache);
    param[4] = 0.5*param[3]/w;
    param[5] = 0.2*w;

    return ok && w >= 0.0;
}

static gboolean
smthpyr4_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                  FitShapeEstimateCache *estimcache)
{
    gdouble w, a_unused;
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_estimate(xyz, n,
                                      param + 0, param + 1, param + 2, param + 3, &w, &a_unused, param + 6,
                                      estimcache);
    param[4] = 0.5*param[3]/w;
    param[5] = 0.2*w;

    return ok && w >= 0.0;
}

/*******************************************************************************************************************
 *
 * Gaussian broadened 3-sided pyramid
 *
 *******************************************************************************************************************/

/* Approximation of the difference between sharp and smooth pyramid, in
 * suitable coordinates transformed to [0,1]². */
static inline double
ratapprox_p3(double u, double v)
{
    double a = 1.03624110601701;
    double B = -0.024290223394786;
    double b = -0.364959025944882;
    double c = -0.374337838348571;
    double C = 0.407619604407758;
    double d = 0.289240769832989;
    double D = 0.640668456580655;
    double e = 0.312954075690613;
    double E = -0.313229623765805;
    double F = 0.0384356719145661;
    double f = -0.176113198624078;
    double g = 0.0121632908390103;
    double G = -0.129217472902444;
    double h = 0.254015515724361;
    double H = 0.393039175617963;
    double v2 = v*v;
    double u2 = u*u;
    double u4 = u2*u2;
    double numer = (1 + b*u + c*v + d*u*v + e*v2 + f*u2 + g*u4 + h*u*v2 + G*u2*v2);
    double denom = (1 + C*v + D*u*v + E*u2 + H*u*v2 + B*u4 + F*u2*v2);

    return a*(1 - u)*numer/denom;
}

static inline double
infapprox_p3(double x, double y)
{
    double a = 1.03624110601701;
    double u, v, xx;

    /* move origin of x to the line y = x√3 */
    xx = x - y/GWY_SQRT3;

    /* log-transform from [0,∞) to [0,1) */
    u = sqrt(1 - exp(-xx*xx));
    v = sqrt(1 - exp(-y*y));

    return ratapprox_p3(u, v) - a + M_SQRT2*x;
}

static gdouble
smthpyr3_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble a = fabs(param[4]);
    gdouble R = fabs(param[5]);
    gdouble phi = param[6];
    gdouble bx = param[7];
    gdouble by = param[8];
    gdouble sigmaq2, val, t, cphi, sphi;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    if (G_UNLIKELY(!h || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    /* fold coordinates to the unique sextant of the plane */
    y = fabs(y);
    if (y > x*GWY_SQRT3) {
        t = (-x + y*GWY_SQRT3)/2;
        y = fabs((-x*GWY_SQRT3 - y)/2);
        x = t;
    }

    if (G_UNLIKELY(!R)) {
        val = x*a;
    }
    else {
        sigmaq2 = 0.75*GWY_SQRT3/GWY_SQRT_PI*R*a;
        val = infapprox_p3(x/sigmaq2, y/sigmaq2);
        val *= sigmaq2*a/G_SQRT2;
    }

    if (val >= fabs(h))
        return z0;

    return h > 0.0 ? z0 + (h - val) : z0 - (h - val);
}

static gboolean
smthpyr3_init(const GwyXYZ *xyz, guint n, gdouble *param,
              FitShapeEstimateCache *estimcache)
{
    gdouble w, a_unused;
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_init(xyz, n,
                                  param + 0, param + 1, param + 2, param + 3, &w, &a_unused, param + 6,
                                  estimcache);
    param[4] = 0.5*param[3]/w;
    param[5] = 0.2*w;

    return ok && w >= 0.0;
}

static gboolean
smthpyr3_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                  FitShapeEstimateCache *estimcache)
{
    gdouble w, a_unused;
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    ok = common_bump_feature_estimate(xyz, n,
                                      param + 0, param + 1, param + 2, param + 3, &w, &a_unused, param + 6,
                                      estimcache);
    param[4] = 0.5*param[3]/w;
    param[5] = 0.2*w;

    return ok && w >= 0.0;
}

/*******************************************************************************************************************
 *
 * Cylinder
 *
 *******************************************************************************************************************/

static gdouble
cylinder_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble h = param[3];
    gdouble R = fabs(param[4]);
    gdouble C = fabs(param[5]);
    gdouble a = fabs(param[6]);
    gdouble phi = param[7];
    gdouble bx = param[8];
    gdouble by = param[9];
    gdouble t, val, cphi, sphi;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    if (G_UNLIKELY(!R || !a))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t;

    t = (x*x*a + y*y/a)/(R*R);
    if (t == 0.0)
        val = z0 + h;
    else if (C == 0.0) {
        val = z0 + (t < 1.0 ? h : 0.0);
    }
    else {
        t = sqrt(t);
        t = t - 1.0/t;
        val = z0 + 0.5*h*(1.0 - erf(t/C));
    }

    return val;
}

static void
cylinder_diff(gdouble x, gdouble y,
              const gdouble *param,
              const gboolean *fixed_param,
              gdouble *der)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    /*gdouble z0 = param[2];*/
    gdouble h = param[3];
    gdouble R = fabs(param[4]);
    gdouble C = fabs(param[5]);
    gdouble a = fabs(param[6]);
    gdouble phi = param[7];
    gdouble bx = param[8];
    gdouble by = param[9];
    gdouble t2, u, xx, yy, cphi, sphi, E, dudt, aa;

    if (G_UNLIKELY(a == 0.0 || R == 0.0 || C == 0.0)) {
        gwy_clear(der, G_N_ELEMENTS(cylinder_params));
        return;
    }

    x -= xc;
    y -= yc;

    der[2] = 1;
    der[8] = x;
    der[9] = y;

    _gwy_sincos(phi, &sphi, &cphi);
    xx = x*cphi - y*sphi;
    yy = x*sphi + y*cphi;
    t2 = xx*xx*a + yy*yy/a;

    if (G_UNLIKELY(t2 == 0.0)) {
        der[0] = der[1] = der[4] = der[5] = der[6] = der[7] = 0.0;
        der[3] = 1.0;
        return;
    }

    u = sqrt(t2);
    dudt = 0.5/u*(1.0/R + R/t2);
    u = (u/R - R/u)/C;
    E = exp(-u*u);
    aa = a - 1.0/a;

    if (!fixed_param || !fixed_param[3])
        der[3] = 0.5*(1.0 - erf(u));

    der[0] = -bx + 2.0*h/GWY_SQRT_PI*E/C * dudt * (a*xx*cphi + yy/a*sphi);
    der[1] = -by - 2.0*h/GWY_SQRT_PI*E/C * dudt * (a*xx*sphi - yy/a*cphi);
    der[4] = h/GWY_SQRT_PI*E/C * (sqrt(t2)/(R*R) + 1.0/sqrt(t2));
    der[5] = h/GWY_SQRT_PI*E/C * u;
    der[6] = h/GWY_SQRT_PI*E/C * dudt * (yy*yy/(a*a) - xx*xx);
    der[7] = 2.0*h/GWY_SQRT_PI*E/C * dudt * aa*xx*yy;
    /*
    compare_der(x + xc, y + yc, param, der,
                G_N_ELEMENTS(cylinder_params), cylinder_func,
                (1<<2)|(1<<3)|(1<<4)|(1<<5)|(1<<6)|(1<<7)|(1<<8)|(1<<9));
                */
}

static gboolean
cylinder_init(const GwyXYZ *xyz, guint n, gdouble *param,
              FitShapeEstimateCache *estimcache)
{
    gboolean ok;

    param[8] = 0.0;
    param[9] = 0.0;
    param[5] = 0.05;
    ok = common_bump_feature_init(xyz, n,
                                  param + 0, param + 1, param + 2, param + 3, param + 4, param + 6, param + 7,
                                  estimcache);

    return ok && param[4] >= 0.0;
}

static gboolean
cylinder_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                  FitShapeEstimateCache *estimcache)
{
    gboolean ok;

    param[7] = 0.0;
    param[8] = 0.0;
    param[5] = 0.05;
    ok = common_bump_feature_estimate(xyz, n,
                                      param + 0, param + 1, param + 2, param + 3, param + 4, param + 6, param + 7,
                                      estimcache);

    return ok && param[4] >= 0.0;
}

static gdouble
cylinder_calc_R1(const gdouble *param)
{
    return param[4]/sqrt(fabs(param[6]));
}

static gdouble
cylinder_calc_err_R1(const gdouble *param,
                     const gdouble *param_err,
                     const gdouble *correl)
{
    gdouble R1 = cylinder_calc_R1(param);
    gdouble diff[G_N_ELEMENTS(cylinder_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[4] = R1/param[4];
    diff[6] = -0.5*R1/param[6];
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

static gdouble
cylinder_calc_R2(const gdouble *param)
{
    return param[4]*sqrt(fabs(param[6]));
}

static gdouble
cylinder_calc_err_R2(const gdouble *param,
                     const gdouble *param_err,
                     const gdouble *correl)
{
    gdouble R1 = cylinder_calc_R1(param);
    gdouble diff[G_N_ELEMENTS(cylinder_params)];

    gwy_clear(diff, G_N_ELEMENTS(diff));
    diff[4] = R1/param[4];
    diff[6] = 0.5*R1/param[6];
    return dotprod_with_correl(diff, param_err, correl, G_N_ELEMENTS(diff));
}

/*******************************************************************************************************************
 *
 * Lunette
 *
 *******************************************************************************************************************/

static gdouble
lunette_func(gdouble x, gdouble y, const gdouble *param)
{
    gdouble xc = param[0];
    gdouble yc = param[1];
    gdouble z0 = param[2];
    gdouble R = fabs(param[3]);
    gdouble phi = param[4];
    gdouble h = fabs(param[5]);
    gdouble delta = fabs(param[6]);
    gdouble r = fabs(param[7]);
    gdouble bx = param[8];
    gdouble by = param[9];
    gdouble p = fabs(param[10]);
    gdouble t1, t2, val, cphi, sphi, m, R2 = R*R;

    x -= xc;
    y -= yc;
    z0 += bx*x + by*y;

    if (G_UNLIKELY(!R))
        return z0;

    _gwy_sincos(phi, &sphi, &cphi);
    t1 = x*cphi - y*sphi;
    y = x*sphi + y*cphi;
    x = t1;

    t1 = x*x + y*y;
    t2 = (x + delta)*(x + delta) + y*y;
    val = z0;
    if (t2 < R2) {
        val += h*p*sqrt(R2 - t2)/R;
        if (t1 > r*r) {
            m = fmin(fabs(R + delta - r), 2*R);
            t1 = r - sqrt(t1);
            t2 = sqrt(t2) - R;
            val += h*t1*t2/(0.25*m*m);
        }
    }

    return val;
}

static gboolean
lunette_init(const GwyXYZ *xyz, guint n, gdouble *param,
             FitShapeEstimateCache *estimcache)
{
    gdouble xm, ym, r, zmin, zmax;

    circumscribe_x_y(xyz, n, &xm, &ym, &r, estimcache);
    range_z(xyz, n, &zmin, &zmax, estimcache);

    param[0] = xm;
    param[1] = ym;
    param[2] = zmin;
    param[3] = r/2.0;
    param[4] = 0.0;
    param[5] = zmax - zmin;
    param[6] = r/15.0;
    param[7] = r/3.0;
    param[8] = 0.0;
    param[9] = 0.0;
    param[10] = 0.2;

    return TRUE;
}

static gboolean
lunette_estimate(const GwyXYZ *xyz, guint n, gdouble *param,
                 FitShapeEstimateCache *estimcache)
{
    /* Note that 16 bins correspond to 8 directions.
     * Some bit combinations cannot occur, these are marked -1. */
    enum { NBINS = 16 };
    static gint phitable[NBINS] = { 4, -1, 3, -1, 5, 6, -1, -1, -1, -1, 2, 1, -1, 7, -1, 0 };
    gdouble zsum[NBINS];
    guint wsum[NBINS];
    gdouble xm, ym, r;
    guint i, imax;

    lunette_init(xyz, n, param, estimcache);
    circumscribe_x_y(xyz, n, &xm, &ym, &r, estimcache);
    gwy_clear(zsum, NBINS);
    gwy_clear(wsum, NBINS);
    for (i = 0; i < n; i++) {
        gdouble x = xyz[i].x - xm, y = xyz[i].y - ym;
        guint idx = ((x < 0.0) | ((y > 0.0) << 1)
                     | ((-x > y) << 2) | ((x < y) << 3));
        zsum[idx] += xyz[i].z;
        wsum[idx]++;
    }
    imax = 0;
    for (i = 0; i < NBINS; i++) {
        gwy_debug("[%u] %g %d (%d)", i, zsum[i], wsum[i], phitable[i]);
        if (!wsum[i])
            continue;
        if (!wsum[imax] || zsum[i]/wsum[i] > zsum[imax]/wsum[imax])
            imax = i;
    }
    g_assert(phitable[imax] != -1);
    param[4] = G_PI/4.0*(phitable[imax] + 0.5);
    gwy_debug("imax=%u, phi = %g deg", imax, 180.0/G_PI*param[4]);

    return TRUE;
}

/*******************************************************************************************************************
 *
 * Documentation
 *
 *******************************************************************************************************************/

/**
 * SECTION:gwyshapefitpreset
 * @title: GwyShapeFitPreset
 * @short_description: 3D geometrical shape fitting functions
 * @see_also: #GwyNLFitter
 *
 * <link linkend="GwyNLFitter">Non-linear fitter</link> 3D geometrical shape presets are predefined fitting functions
 * that can be used to determine the parameters of various features in images and XYZ data.
 *
 * As of version 2.47 the defined functions include:
 * <simplelist type='vert'>
 * <member><literal>"Grating (simple)"</literal></member>
 * <member><literal>"Grating (3-level)"</literal></member>
 * <member><literal>"Holes"</literal></member>
 * <member><literal>"Sphere"</literal></member>
 * <member><literal>"Cylinder (lying)"</literal></member>
 * <member><literal>"Gaussian"</literal></member>
 * <member><literal>"Lorentzian"</literal></member>
 * <member><literal>"Pyramid (diamond)"</literal></member>
 * <member><literal>"Parabolic bump"</literal></member>
 * <member><literal>"Ring"</literal></member>
 * </simplelist>
 *
 * Since version 2.53 they also include:
 * <simplelist type='vert'>
 * <member><literal>"Cone"</literal></member>
 * <member><literal>"Plane"</literal></member>
 * <member><literal>"Pyramid (3-sided)"</literal></member>
 * <member><literal>"Step (one-sided)"</literal></member>
 * <member><literal>"Step (two-sided)"</literal></member>
 * </simplelist>
 *
 * Since version 2.54 they also include:
 * <simplelist type='vert'>
 * <member><literal>"Cylinder"</literal></member>
 * </simplelist>
 *
 * Since version 2.56 they also include:
 * <simplelist type='vert'>
 * <member><literal>"Half-sphere"</literal></member>
 * </simplelist>
 *
 * Since version 2.57 they also include:
 * <simplelist type='vert'>
 * <member><literal>"Smooth cone"</literal></member>
 * <member><literal>"Smooth pyramid"</literal></member>
 * <member><literal>"Smooth pyramid (3-sided)"</literal></member>
 * </simplelist>
 *
 * Since version 2.59 they also include:
 * <simplelist type='vert'>
 * <member><literal>"Exponential"</literal></member>
 * <member><literal>"Lunette"</literal></member>
 * </simplelist>
 *
 * Since version 2.61 they also include:
 * <simplelist type='vert'>
 * <member><literal>"Pyramid (rectangle)"</literal></member>
 * </simplelist>
 **/

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
