/*
 *  $Id: gwyprocessenums.c 25400 2023-06-02 13:50:51Z yeti-dn $
 *  Copyright (C) 2005-2021 David Necas (Yeti), Petr Klapetek.
 *  E-mail: yeti@gwyddion.net, klapetek@gwyddion.net.
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
#include <libgwyddion/gwymacros.h>
#include <libprocess/gwyprocessenums.h>

/**
 * GwyMergeType:
 * @GWY_MERGE_UNION: Union (logical or) merging.
 * @GWY_MERGE_INTERSECTION: Intersection (logical and) merging.
 *
 * Mask merge type (namely used in grain processing).
 **/
/**
 * gwy_merge_type_get_enum:
 *
 * Returns #GwyEnum for #GwyMergeType enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 **/
const GwyEnum*
gwy_merge_type_get_enum(void)
{
    static const GwyEnum entries[] = {
        { N_("Union"),            GWY_MERGE_UNION,        },
        { N_("Intersection"),     GWY_MERGE_INTERSECTION, },
        { NULL,                   0,                      },
    };
    return entries;
}

/**
 * GwyPlaneSymmetry:
 * @GWY_SYMMETRY_AUTO: Automatic symmetry selection.
 * @GWY_SYMMETRY_PARALLEL: Parallel symmetry, there is one prevalent direction (bilateral).
 * @GWY_SYMMETRY_TRIANGULAR: Triangular symmetry, there are three prevalent directions (unilateral) by 120 degrees.
 * @GWY_SYMMETRY_SQUARE: Square symmetry, two prevalent directions (bilateral) oriented approximately along image
 *                       sides.
 * @GWY_SYMMETRY_RHOMBIC: Rhombic symmetry, two prevalent directions (bilateral) oriented approximately along
 *                        diagonals.
 * @GWY_SYMMETRY_HEXAGONAL: Hexagonal symmetry, three prevalent directions (bilateral) by 120 degrees.
 * @GWY_SYMMETRY_LAST: The number of symmetries.
 *
 * Plane symmetry types for rotation correction.
 **/
/**
 * gwy_plane_symmetry_get_enum:
 *
 * Returns #GwyEnum for #GwyPlaneSymmetry enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 **/
const GwyEnum*
gwy_plane_symmetry_get_enum(void)
{
    static const GwyEnum entries[] = {
        { N_("Detected"),   GWY_SYMMETRY_AUTO       },
        { N_("Parallel"),   GWY_SYMMETRY_PARALLEL   },
        { N_("Triangular"), GWY_SYMMETRY_TRIANGULAR },
        { N_("symmetry|Square"),     GWY_SYMMETRY_SQUARE     },
        { N_("Rhombic"),    GWY_SYMMETRY_RHOMBIC    },
        { N_("Hexagonal"),  GWY_SYMMETRY_HEXAGONAL  },
        { NULL,             0,                      },
    };
    return entries;
};

/**
 * gwy_2d_cwt_wavelet_type_get_enum:
 *
 * Returns #GwyEnum for #Gwy2DCWTWaveletType enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 **/
const GwyEnum*
gwy_2d_cwt_wavelet_type_get_enum(void)
{
    static const GwyEnum entries[] = {
        { N_("Gaussian"),          GWY_2DCWT_GAUSS      },
        { N_("Hat"),               GWY_2DCWT_HAT        },
        { NULL,                    0,                   },
    };
    return entries;
}

/**
 * GwyOrientation:
 * @GWY_ORIENTATION_HORIZONTAL: Horizontal orientation.
 * @GWY_ORIENTATION_VERTICAL: Vertical orientation.
 *
 * Orientation type.
 **/
/**
 * gwy_orientation_get_enum:
 *
 * Returns #GwyEnum for #GwyOrientation enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 **/
const GwyEnum*
gwy_orientation_get_enum(void)
{
    static const GwyEnum entries[] = {
        { N_("Horizontal"),  GWY_ORIENTATION_HORIZONTAL,  },
        { N_("Vertical"),    GWY_ORIENTATION_VERTICAL,    },
        { NULL,              0,                           },
    };
    return entries;
}

/**
 * gwy_dwt_type_get_enum:
 *
 * Returns #GwyEnum for #GwyDWTType enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 **/
const GwyEnum*
gwy_dwt_type_get_enum(void)
{
    static const GwyEnum entries[] = {
        { N_("Haar"),          GWY_DWT_HAAR    },
        { N_("Daubechies 4"),  GWY_DWT_DAUB4   },
        { N_("Daubechies 6"),  GWY_DWT_DAUB6   },
        { N_("Daubechies 8"),  GWY_DWT_DAUB8   },
        { N_("Daubechies 12"), GWY_DWT_DAUB12  },
        { N_("Daubechies 20"), GWY_DWT_DAUB20  },
        { NULL,                0,              },
     };
    return entries;
}

/**
 * gwy_dwt_denoise_type_get_enum:
 *
 * Returns #GwyEnum for #GwyDWTDenoiseType enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 **/
const GwyEnum*
gwy_dwt_denoise_type_get_enum(void)
{
    static const GwyEnum entries[] = {
        { N_("Universal"),                 GWY_DWT_DENOISE_UNIVERSAL,       },
        { N_("Scale adaptive"),            GWY_DWT_DENOISE_SCALE_ADAPTIVE,  },
        { N_("Scale and space adaptive"),  GWY_DWT_DENOISE_SPACE_ADAPTIVE,  },
        { NULL,                            0,                               },
    };
    return entries;
}

/**
 * GwyInterpolationType:
 * @GWY_INTERPOLATION_NONE: No interpolation at all, resulting values are not
 *                          defined, it must not be used for interpolation.
 *                          It can be used in resize operations discarding
 *                          original data.
 * @GWY_INTERPOLATION_ROUND: Round interpolation (more precisely symmetric
 *                           nearest neighbour interpolation).
 * @GWY_INTERPOLATION_LINEAR: Linear interpolation.
 * @GWY_INTERPOLATION_BILINEAR: Old name for %GWY_INTERPOLATION_LINEAR.  This
 *                              is a misnomer because it in fact denotes linear
 *                              interpolation of arbitrary dimension.  Use
 *                              %GWY_INTERPOLATION_LINEAR instead.
 * @GWY_INTERPOLATION_KEY: Cubic Key's interpolation (with a=-1/2).
 * @GWY_INTERPOLATION_BSPLINE: B-spline interpolation.
 * @GWY_INTERPOLATION_OMOMS: Omoms interpolation.
 * @GWY_INTERPOLATION_NNA: Nearest neighbour approximation.
 * @GWY_INTERPOLATION_SCHAUM: Cubic Schaum interpolation.
 *
 * Interpolation types.
 **/
/**
 * gwy_interpolation_type_get_enum:
 *
 * Returns #GwyEnum for #GwyInterpolationType enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 **/
const GwyEnum*
gwy_interpolation_type_get_enum(void)
{
    static const GwyEnum entries[] = {
        /* TRANSLATORS: Interpolation type (AKA nearest neighbour) */
        { N_("Round"),   GWY_INTERPOLATION_ROUND,   },
        { N_("Linear"),  GWY_INTERPOLATION_LINEAR,  },
        { N_("Key"),     GWY_INTERPOLATION_KEY,     },
        { N_("BSpline"), GWY_INTERPOLATION_BSPLINE, },
        { N_("OMOMS"),   GWY_INTERPOLATION_OMOMS,   },
        { N_("NNA"),     GWY_INTERPOLATION_NNA,     },
        { N_("Schaum"),  GWY_INTERPOLATION_SCHAUM,  },
        { NULL,          0,                         },
    };
    return entries;
}

/**
 * GwyWindowingType:
 * @GWY_WINDOWING_NONE: No windowing is applied.
 * @GWY_WINDOWING_HANN: Hann window.
 * @GWY_WINDOWING_HAMMING: Hamming window.
 * @GWY_WINDOWING_BLACKMANN: Blackmann window.
 * @GWY_WINDOWING_LANCZOS: Lanczos window.
 * @GWY_WINDOWING_WELCH: Welch window.
 * @GWY_WINDOWING_RECT: Rectangular window.
 * @GWY_WINDOWING_NUTTALL: Nuttall window (Since 2.7).
 * @GWY_WINDOWING_FLAT_TOP: Flat-top window (Since 2.7).
 * @GWY_WINDOWING_KAISER25: Kaiser window with β=2.5 (Since 2.7).
 *
 * Frequency windowing type.
 **/
/**
 * gwy_windowing_type_get_enum:
 *
 * Returns #GwyEnum for #GwyWindowingType enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 **/
const GwyEnum*
gwy_windowing_type_get_enum(void)
{
    static const GwyEnum entries[] = {
        { N_("windowing|None"), GWY_WINDOWING_NONE,      },
        { N_("Hann"),           GWY_WINDOWING_HANN,      },
        { N_("Hamming"),        GWY_WINDOWING_HAMMING,   },
        { N_("Blackmann"),      GWY_WINDOWING_BLACKMANN, },
        { N_("Lanczos"),        GWY_WINDOWING_LANCZOS,   },
        { N_("Welch"),          GWY_WINDOWING_WELCH,     },
        { N_("Rect"),           GWY_WINDOWING_RECT,      },
        { N_("Nuttall"),        GWY_WINDOWING_NUTTALL,   },
        { N_("Flat-top"),       GWY_WINDOWING_FLAT_TOP,  },
        { N_("Kaiser 2.5"),     GWY_WINDOWING_KAISER25,  },
        { NULL,                 0,                       },
    };
    return entries;
}

/**
 * gwy_correlation_type_get_enum:
 *
 * Returns #GwyEnum for #GwyCorrelationType enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 **/
const GwyEnum*
gwy_correlation_type_get_enum(void)
{
    static const GwyEnum entries[] = {
        { N_("correlation|Normal"), GWY_CORRELATION_NORMAL, },
        { N_("FFT"),                GWY_CORRELATION_FFT,    },
        { N_("Phase only"),         GWY_CORRELATION_POC,    },
        { NULL,                     0,                      },
    };
    return entries;
}

/**
 * GwyMaskingType:
 * @GWY_MASK_EXCLUDE: Exclude data under mask, i.e. take into account only
 *                    data not covered by the mask.
 * @GWY_MASK_INCLUDE: Take into account only data under the mask.
 * @GWY_MASK_IGNORE: Ignore mask, if present, and use all data.
 *
 * Mask handling in procedures that can apply masking.
 *
 * Note at present many procedures do not have a masking argument and hence cannot apply masks in exclude mode.
 *
 * Since: 2.12
 **/
/**
 * gwy_masking_type_get_enum:
 *
 * Returns #GwyEnum for #GwyMaskingType enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 *
 * Since: 2.18
 **/
const GwyEnum*
gwy_masking_type_get_enum(void)
{
    static const GwyEnum entries[] = {
        { N_("Exclude masked region"),          GWY_MASK_EXCLUDE, },
        { N_("Include only masked region"),     GWY_MASK_INCLUDE, },
        { N_("Use entire image (ignore mask)"), GWY_MASK_IGNORE,  },
        { NULL,                                 0,                },
    };
    return entries;
}


/************************** Documentation ****************************/

/**
 * SECTION:gwyprocessenums
 * @title: gwyprocessenums
 * @short_description: Common enumerations
 **/

/**
 * GwyTransformDirection:
 * @GWY_TRANSFORM_DIRECTION_BACKWARD: Backward (inverse) transform.
 * @GWY_TRANSFORM_DIRECTION_FORWARD: Forward (direct) transform.
 *
 * Transform (namely integral transform) direction.
 *
 * In FFT, it is equal to sign of the exponent, that is the backward transform uses -1, the forward transform +1.
 * This is the opposite sign convention to FFTW (for instance), so care must be taken when mixing operations.
 **/

/**
 * GwyPlaneFitQuantity:
 * @GWY_PLANE_FIT_A: Constant coefficient (mean value).
 * @GWY_PLANE_FIT_BX: Linear coefficient in @x, where @x is in pixel coordinates.
 * @GWY_PLANE_FIT_BY: Linear coefficient in @y, where @y is in pixel coordinates.
 * @GWY_PLANE_FIT_ANGLE: Slope orientation in (x,y) plane (in radians).
 * @GWY_PLANE_FIT_SLOPE: Absolute slope value (that is sqrt(bx*bx + by*by)).
 * @GWY_PLANE_FIT_S0: Residual sum of squares.
 * @GWY_PLANE_FIT_S0_REDUCED: Slope-reduced residual sum of squares.
 *
 * Local plane fitting quantity to request from gwy_data_field_area_fit_local_planes() and similar functions.
 **/

/**
 * GwyGrainQuantity:
 * @GWY_GRAIN_VALUE_PROJECTED_AREA: Projected (flat) grain area.
 * @GWY_GRAIN_VALUE_EQUIV_SQUARE_SIDE: Side of a square with the same area as the grain.
 * @GWY_GRAIN_VALUE_EQUIV_DISC_RADIUS: Radius of a disc with the same area as the grain.
 * @GWY_GRAIN_VALUE_SURFACE_AREA: Surface area.
 * @GWY_GRAIN_VALUE_MAXIMUM: Minimum value.
 * @GWY_GRAIN_VALUE_MINIMUM: Maximum value.
 * @GWY_GRAIN_VALUE_MEAN: Mean value.
 * @GWY_GRAIN_VALUE_MEDIAN: Median value.
 * @GWY_GRAIN_VALUE_PIXEL_AREA: Flat grain area measured in pixels.  This value is redundant but it is useful for
 *                              filtering (Since 2.37).
 * @GWY_GRAIN_VALUE_HALF_HEIGHT_AREA: Projected area of the part of grain that is above the half-height, i.e. the
 *                                    height between the minimum and maximum (Since 2.8).
 * @GWY_GRAIN_VALUE_RMS: Standard deviation of grain values.  (Since 2.51)
 * @GWY_GRAIN_VALUE_FLAT_BOUNDARY_LENGTH: Length of projected grain boundary. Note for grains not topologically
 *                                        equivalent to disc, only the length of the outer boundary is calculated.
 * @GWY_GRAIN_VALUE_MINIMUM_BOUND_SIZE: Minimum lateral bounding size, in other words the minimum length of grain
 *                                      projection to any line in the horizontal plane.
 * @GWY_GRAIN_VALUE_MINIMUM_BOUND_ANGLE: Direction of the minimum lateral bounding size (arbitrary one if the minimum
 *                                       is not unique).
 * @GWY_GRAIN_VALUE_MAXIMUM_BOUND_SIZE: Maximum lateral bounding size, in other words the maximum length of grain
 *                                      projection to any line in the horizontal plane.
 * @GWY_GRAIN_VALUE_MAXIMUM_BOUND_ANGLE: Direction of the maximum lateral bounding size (arbitrary one if the maximum
 *                                       is not unique).
 * @GWY_GRAIN_VALUE_CENTER_X: Grain centre horizontal position, i.e. the mean value of its physical x-coordinates.
 *                            (Since 2.7)
 * @GWY_GRAIN_VALUE_CENTER_Y: Grain centre vertical position, i.e. the mean value of its physical y-coordinates.
 *                            (Since 2.7)
 * @GWY_GRAIN_VALUE_VOLUME_0: Grain volume calculated with grain basis at @z=0 (therefore it is just an integral it can
 *                            be negative).  (Since 2.3)
 * @GWY_GRAIN_VALUE_VOLUME_MIN: Grain volume calculated with grain basis at grain minimum value.  This value is
 *                              a lower bound.  (Since 2.3)
 * @GWY_GRAIN_VALUE_VOLUME_LAPLACE: Grain volume calculated with grain basis calculated by laplacian interpolation of
 *                                  surrounding values.  (Since 2.3)
 * @GWY_GRAIN_VALUE_SLOPE_THETA: Spherical angle theta of grain normal (0 is upwards). (Since 2.7)
 * @GWY_GRAIN_VALUE_SLOPE_PHI: Spherical angle phi of grain normal (0 is in positive x direction). (Since 2.7)
 * @GWY_GRAIN_VALUE_BOUNDARY_MINIMUM: Minimum value on the grain inner boundary.  (Since 2.12)
 * @GWY_GRAIN_VALUE_BOUNDARY_MAXIMUM: Maximum value in the grain inner boundary.  (Since 2.12)
 * @GWY_GRAIN_VALUE_CURVATURE_CENTER_X: Grain curvature centre horizontal position.  For too small or flat grains it
 *                                      reduces to the horizontal position of geometrical centre. (Since 2.22)
 * @GWY_GRAIN_VALUE_CURVATURE_CENTER_Y: Grain curvature centre vertical position.  For too small or flat grains it
 *                                      reduces to the vertical position of geometrical centre. (Since 2.22)
 * @GWY_GRAIN_VALUE_CURVATURE_CENTER_Z: The value at curvature centre.  Note this is the value in the origin of the
 *                                      fitted quadratic surface, not at the real surface. (Since 2.22)
 * @GWY_GRAIN_VALUE_CURVATURE1: Smaller grain curvature. (Since 2.22)
 * @GWY_GRAIN_VALUE_CURVATURE2: Larger grain curvature. (Since 2.22)
 * @GWY_GRAIN_VALUE_CURVATURE_ANGLE1: Direction of the smaller grain curvature radius.  If the grain is flat or too
 *                                    small the angle is reported as 0. (Since 2.22)
 * @GWY_GRAIN_VALUE_CURVATURE_ANGLE2: Direction of the larger grain curvature radius.  If the grain is flat or too
 *                                    small the angle is reported as π/2. (Since 2.22)
 * @GWY_GRAIN_VALUE_INSCRIBED_DISC_R: Radius of maximum disc that fits inside the grain (Since 2.30)
 * @GWY_GRAIN_VALUE_INSCRIBED_DISC_X: Real X-coordinate of the centre of the maximum inscribed disc. (Since 2.30)
 * @GWY_GRAIN_VALUE_INSCRIBED_DISC_Y: Real Y-coordinate of the centre of the maximum inscribed disc. (Since 2.30)
 * @GWY_GRAIN_VALUE_CONVEX_HULL_AREA: Projected (flat) area of grain convex hull. (Since 2.30)
 * @GWY_GRAIN_VALUE_CIRCUMCIRCLE_R: Radius of minimum circle containing the grain.  (Since 2.30)
 * @GWY_GRAIN_VALUE_CIRCUMCIRCLE_X: Real X-coordinate of the centre of the minimum circumcircle. (Since 2.30)
 * @GWY_GRAIN_VALUE_CIRCUMCIRCLE_Y: Real Y-coordinate of the centre of the minimum circumcircle. (Since 2.30)
 * @GWY_GRAIN_VALUE_MEAN_RADIUS: Mean distance from boundary to the grain centre as defined by
 *                               @GWY_GRAIN_VALUE_CENTER_X and @GWY_GRAIN_VALUE_CENTER_Y. (Since 2.30)
 * @GWY_GRAIN_VALUE_EQUIV_ELLIPSE_MAJOR: Length of major semiaxis of equivalent ellipse. (Since 2.36)
 * @GWY_GRAIN_VALUE_EQUIV_ELLIPSE_MINOR: Length of minor semiaxis of equivalent ellipse. (Since 2.36)
 * @GWY_GRAIN_VALUE_EQUIV_ELLIPSE_ANGLE: Orientation of the major axis of equivalent ellipse. (Since 2.36)
 * @GWY_GRAIN_VALUE_MINIMUM_MARTIN_DIAMETER: Minimum value of Martin diameter. (Since 2.50)
 * @GWY_GRAIN_VALUE_MINIMUM_MARTIN_ANGLE: Direction corresponding to minimum Martin diameter.  (Since 2.50)
 * @GWY_GRAIN_VALUE_MAXIMUM_MARTIN_DIAMETER: Maximum value of Martin diameter. (Since 2.50)
 * @GWY_GRAIN_VALUE_MAXIMUM_MARTIN_ANGLE: Direction corresponding to maximum Martin diameter.  (Since 2.50)
 *
 * Grain quantity to request from gwy_data_field_grains_get_distribution() and similar functions.
 **/

/**
 * GwyDataCompatibilityFlags:
 * @GWY_DATA_COMPATIBILITY_RES: Pixel sizes.
 * @GWY_DATA_COMPATIBILITY_REAL: Real (physical) dimensions.
 * @GWY_DATA_COMPATIBILITY_MEASURE: Real to pixel ratios.
 * @GWY_DATA_COMPATIBILITY_LATERAL: Units of lateral dimensions.
 * @GWY_DATA_COMPATIBILITY_VALUE: Units of values (for all curves in the case of lawns).
 * @GWY_DATA_COMPATIBILITY_AXISCAL: Axis calibrations.  At present it only makes sense for #GwyBrick which can have
 *                                  Z-calibrations.  (Since 2.51)
 * @GWY_DATA_COMPATIBILITY_NCURVES: Number of lawn curves.  (Since 2.60)
 * @GWY_DATA_COMPATIBILITY_CURVELEN: Lengths of curves in all lawn pixels.  (Since 2.60)
 * @GWY_DATA_COMPATIBILITY_ALL: Mask of all defined flags.
 *
 * Data line, field, brick and lawn compatibility flags.
 *
 * It is not recommended to pass %GWY_DATA_COMPATIBILITY_ALL to checking functions since not all flags are meaningful
 * for all data objects (even though meaningless flags are generally silently ignored).
 **/

/**
 * GwyDataFieldCached:
 * @GWY_DATA_FIELD_CACHE_MIN: Overall minimum.
 * @GWY_DATA_FIELD_CACHE_MAX: Overall maximum.
 * @GWY_DATA_FIELD_CACHE_SUM: Sum of all values.
 * @GWY_DATA_FIELD_CACHE_RMS: Root mean square.
 * @GWY_DATA_FIELD_CACHE_MED: Median.
 * @GWY_DATA_FIELD_CACHE_ARF: Auto-range from.
 * @GWY_DATA_FIELD_CACHE_ART: Auto-range to.
 * @GWY_DATA_FIELD_CACHE_ARE: Surface area.
 * @GWY_DATA_FIELD_CACHE_VAR: Variation.
 * @GWY_DATA_FIELD_CACHE_ENT: Entropy.
 * @GWY_DATA_FIELD_CACHE_MSQ: Mean square.
 * @GWY_DATA_FIELD_CACHE_SIZE: The size of statistics cache.
 *
 * Cached data field quantity type.
 *
 * There should be little need to this enum directly except in libgwyprocess
 * methods.
 **/

/**
 * GwyLineStatQuantity:
 * @GWY_LINE_STAT_MEAN: Mean value.
 * @GWY_LINE_STAT_MEDIAN: Median.
 * @GWY_LINE_STAT_MINIMUM: Minimum value.
 * @GWY_LINE_STAT_MAXIMUM: Maximum value.
 * @GWY_LINE_STAT_RMS: Root mean square of deviations from the mean value.
 * @GWY_LINE_STAT_LENGTH: Line length.
 * @GWY_LINE_STAT_SLOPE: Overall line slope.
 * @GWY_LINE_STAT_TAN_BETA0: Root mean square slope.
 * @GWY_LINE_STAT_RA: Arithmetic mean surface roughness.
 * @GWY_LINE_STAT_RZ: Maximum height of the roughness profile.
 * @GWY_LINE_STAT_RT: Total height of the roughness profile.
 * @GWY_LINE_STAT_SKEW: Line skew.
 * @GWY_LINE_STAT_KURTOSIS: Line excess kurtosis (which is 0 for a Gaussaian distrubution, not 3).
 * @GWY_LINE_STAT_RANGE: Difference between maximum and minimum value (Since 2.42).
 * @GWY_LINE_STAT_VARIATION: Variation (integral of absolute value) (Since 2.42).
 * @GWY_LINE_STAT_MINPOS: Minimum position along the line (Since 2.48).
 * @GWY_LINE_STAT_MAXPOS: Maximum position along the line (Since 2.48).
 *
 * Line statistical quantities to be requested with gwy_data_field_area_get_line_stats().
 *
 * Since: 2.2
 **/

/**
 * GwyExteriorType:
 * @GWY_EXTERIOR_UNDEFINED: The values corresponding to or calculated from exterior data values are undefined, they
 *                          may be left unset or set to bogus values.  The caller must handle them itself afterwards,
 *                          for instance by resizing the result to consist of valid data only.
 * @GWY_EXTERIOR_BORDER_EXTEND: Values of exterior pixels are considered to be equal to the values of the nearest
 *                              interior pixels.
 * @GWY_EXTERIOR_MIRROR_EXTEND: The data is considered to be periodically repeated, with odd instances reflected (the
 *                              total period is thus twice the size of the data).
 * @GWY_EXTERIOR_PERIODIC: The data is considered to be periodically repeated.
 * @GWY_EXTERIOR_FIXED_VALUE: Values of exterior pixels are considered to be all equal to a user-specified value.
 * @GWY_EXTERIOR_LAPLACE: Values of exterior pixels are extended using Laplace interpolation (actually extrapolation)
 *                        like gwy_data_field_laplace_solve().  Only some functions implement this method.  (Since
 *                        2.50)
 *
 * Methods to handle pixels outside data.
 *
 * Many methods currently use a fixed metod of handling of exterior pixels, for example area calculation uses
 * extension (border and mirror coincide), convolution uses mirror extension, rotation fills exterior with a fixed
 * value.
 *
 * Since: 2.2
 **/

/**
 * GwyDistanceTransformType:
 * @GWY_DISTANCE_TRANSFORM_CITYBLOCK: City-block distance (sum of horizontal and vertical distances).
 * @GWY_DISTANCE_TRANSFORM_CONN4: Four-connectivity distance; another name for city-block distance.
 * @GWY_DISTANCE_TRANSFORM_CHESS: Chessboard distance (maximum of horizontal and vertical distance).
 * @GWY_DISTANCE_TRANSFORM_CONN8: Eight-connectivity distance; another name for chessboard distance.
 * @GWY_DISTANCE_TRANSFORM_OCTAGONAL48: Octagonal distance beginning from city-block.
 * @GWY_DISTANCE_TRANSFORM_OCTAGONAL84: Octagonal distance beginning from chess.
 * @GWY_DISTANCE_TRANSFORM_OCTAGONAL: Average octagonal distance, i.e. the mean of the 48 and 84 distances (Since
 *                                    2.43).
 * @GWY_DISTANCE_TRANSFORM_EUCLIDEAN: True Euclidean distance (Since 2.43).
 *
 * Type of distance transform.
 *
 * Since: 2.41
 **/
/**
 * gwy_distance_transform_type_get_enum:
 *
 * Returns #GwyEnum for #GwyDistanceTransformType enum type.
 *
 * Returns: %NULL-terminated #GwyEnum which must not be modified nor freed.
 *
 * Since: 2.43
 **/
const GwyEnum*
gwy_distance_transform_type_get_enum(void)
{
    static const GwyEnum entries[] = {
        { N_("distance|City-block"),    GWY_DISTANCE_TRANSFORM_CITYBLOCK,   },
        { N_("distance|Chess"),         GWY_DISTANCE_TRANSFORM_CHESS,       },
        { N_("distance|Octagonal 4,8"), GWY_DISTANCE_TRANSFORM_OCTAGONAL48, },
        { N_("distance|Octagonal 8,4"), GWY_DISTANCE_TRANSFORM_OCTAGONAL84, },
        { N_("distance|Octagonal"),     GWY_DISTANCE_TRANSFORM_OCTAGONAL,   },
        { N_("distance|Euclidean"),     GWY_DISTANCE_TRANSFORM_EUCLIDEAN,   },
        { NULL,                         0,                                  },
    };
    return entries;
}

/**
 * GwyTipType:
 * @GWY_TIP_PYRAMID: N-sided pyramidal tip.
 * @GWY_TIP_PYRAMIDE: Legacy name for %GWY_TIP_PYRAMID.
 * @GWY_TIP_CONTACT: Four-sided pyramidal tip.
 * @GWY_TIP_NONCONTACT: Three-sided pyramidal tip.
 * @GWY_TIP_DELTA: Delta function (single-pixel tip).
 * @GWY_TIP_PARABOLA: Parabolic tip.  (Since 2.45)
 * @GWY_TIP_CONE: Conical tip.  (Since 2.45)
 * @GWY_TIP_ELLPARABOLA: Elliptical parabola tip.  (Since 2.47)
 * @GWY_TIP_BALL_ON_STICK: Spherical tip at the end of almost-cyliner. (Since 2.50)
 *
 * Type of tip shape presets.
 **/

/**
 * GwyTipParamType:
 * @GWY_TIP_PARAM_HEIGHT: Total tip height.  This is used only in the delta function tip; for all others it is
 *                        implied.
 * @GWY_TIP_PARAM_RADIUS: Radius of curvature of the tip apex.
 * @GWY_TIP_PARAM_ROTATION: Rotation angle.
 * @GWY_TIP_PARAM_NSIDES: Number of sides for pyramidal tips.
 * @GWY_TIP_PARAM_SLOPE: Half-angle of the apex (complement of the side slope for straight sides).
 * @GWY_TIP_PARAM_ANISOTROPY: Ratio between larger and smaller tip width in two orthotonal directions.
 *
 * Type of tip model preset parameter.
 *
 * This enum is used with the new tip preset functions gwy_tip_model_get_preset_params(),
 * gwy_tip_model_preset_create(), gwy_tip_model_preset_create_for_zrange().
 *
 * Since: 2.47
 **/

/**
 * GwyMinMaxFilterType:
 * @GWY_MIN_MAX_FILTER_MINIMUM: Minimum filter, i.e. minimum of the surrounding values.
 * @GWY_MIN_MAX_FILTER_EROSION: Another name for the minimum filter.
 * @GWY_MIN_MAX_FILTER_MAXIMUM: Maximum filter, i.e. maximum of the surrounding values.
 * @GWY_MIN_MAX_FILTER_DILATION: Another name for the maximum filter.
 * @GWY_MIN_MAX_FILTER_OPENING: Morphological opening filter.
 * @GWY_MIN_MAX_FILTER_CLOSING: Morphological closing filter.
 * @GWY_MIN_MAX_FILTER_RANGE: Difference between maximum and minimum.
 * @GWY_MIN_MAX_FILTER_NORMALIZATION: Data value rescaled to the range between minimum and maximum.
 *
 * Type of operation based on morphological filters with flat structuring
 * elements.
 *
 * Since: 2.43
 **/

/**
 * GwyRotateResizeType:
 * @GWY_ROTATE_RESIZE_SAME_SIZE: The result has the same area as the original.
 * @GWY_ROTATE_RESIZE_EXPAND: The result is sufficiently large so that all original data are present (usually meaning
 *                            that there will be also lot of exterior).
 * @GWY_ROTATE_RESIZE_CUT: The result is optimally cut to interior data only.
 *
 * Type of rotated data field size determination method.
 *
 * Since: 2.46
 **/

/**
 * GwyAffineScalingType:
 * @GWY_AFFINE_SCALING_AS_GIVEN: Correct lattice vectors lengths are taken as given.
 * @GWY_AFFINE_SCALING_PRESERVE_AREA: Correct lattice vectors are scaled to make the transformation area-preserving.
 * @GWY_AFFINE_SCALING_PRESERVE_X: Correct lattice vectors are scaled to preserve the scale along @x-axis.
 *
 * Type of lattice vector scaling in affine transform preparation.
 *
 * Since: 2.49
 **/

/**
 * GwyCorrSearchType:
 * @GWY_CORR_SEARCH_COVARIANCE_RAW: Raw average of data values multiplied by kernel values.
 * @GWY_CORR_SEARCH_COVARIANCE: Local mean value is subtracted from data before kernel multiplication.
 * @GWY_CORR_SEARCH_COVARIANCE_SCORE: In addition, result is normalised by dividing by the local variance.
 * @GWY_CORR_SEARCH_HEIGHT_DIFF_RAW: Raw mean square difference between data and kernel values.
 * @GWY_CORR_SEARCH_HEIGHT_DIFF: Mean values of data and kernel are adjusted before summing the squared differences.
 * @GWY_CORR_SEARCH_HEIGHT_DIFF_SCORE: In addition, result is normalised by dividing by the local variance.
 * @GWY_CORR_SEARCH_PHASE_ONLY_SCORE: Score from phase-only correlation (since 2.63).
 *
 * Type of correlation search output.
 *
 * Since: 2.50
 **/

/**
 * GwyBrickTransposeType:
 * @GWY_BRICK_TRANSPOSE_XYZ: No change (useful with axis flipping).
 * @GWY_BRICK_TRANSPOSE_XZY: Axes Z and Y are swapped.
 * @GWY_BRICK_TRANSPOSE_YXZ: Axes Y and Z are swapped.
 * @GWY_BRICK_TRANSPOSE_YZX: Axis X becomes Y, Y becomes Z and Z becomes X.
 * @GWY_BRICK_TRANSPOSE_ZXY: Axes X and Z are swapped.
 * @GWY_BRICK_TRANSPOSE_ZYX: Axis X becomes Z, Y becomes X and Z becomes Y.
 *
 * Type of volume data transposition.
 *
 * The enum values names spell which old axis becomes which new axes.
 *
 * Since: 2.51
 **/

/**
 * GwyComputationStateType:
 * @GWY_COMPUTATION_STATE_INIT: Iterator was set up, the next step will actually create temporary data structures and
 *                              precalculate values.
 * @GWY_COMPUTATION_STATE_ITERATE: Iteration is in progress, the @fraction field of state struct holds the fraction
 *                                 completed.
 * @GWY_COMPUTATION_STATE_FINISHED: Calculation has been finished, further calls to the iterator will be no-op.
 *
 * Common iterative computation iterator state type.
 **/

/**
 * GwyComputationState:
 * @state: Current computation state, usually of #GwyComputationStateType, but particular iterators can define their
 *         own types.
 * @fraction: Fraction of computation completed.  For staged algorithms, the fraction refers to the current stage
 *            only.
 *
 * State of iterative computation.
 *
 * Iterators usually append their own private state data, therefore it must not be assumed the public fields @state
 * and @fraction are the only fields.
 *
 * A typical iteration, assuming an iterative computation `foo' with the default #GwyComputationStateType state could
 * be:
 * <informalexample><programlisting>
 * GwyComputationStateType *state;
 * <!-- Hello, gtk-doc! -->
 * state = gwy_data_field_foo_init(GwyDataField *data_field, ...);
 * do {
 *     gwy_data_field_foo_iteration(state);
 *     /<!-- -->* Update progress using state->fraction,
 *         let Gtk+ main loop run, check for abort, ... *<!-- -->/
 *     if (aborted) {
 *         gwy_data_field_foo_finalize(state);
 *         return FALSE;
 *     }
 * } while (state->state != GWY_COMPUTATION_STATE_FINISHED);
 * gwy_data_field_foo_finalize(state);
 * return TRUE;
 * </programlisting></informalexample>
 **/

GType
gwy_computation_state_get_type(void)
{
    static GType compstate_type = 0;

    if (G_UNLIKELY(!compstate_type))
        compstate_type = g_pointer_type_register_static("GwyComputationState");

    return compstate_type;
}

/**
 * gwy_computation_state_get_state:
 * @compstate: Computation state.
 *
 * Gets the state field of a computation state struct.
 *
 * Useful mostly for language bindings.
 *
 * Returns: The state field of @compstate.
 *
 * Since: 2.49
 **/
gint
gwy_computation_state_get_state(GwyComputationState *compstate)
{
    return compstate->state;
}

/**
 * gwy_computation_state_get_fraction:
 * @compstate: Computation state.
 *
 * Gets the fraction field of a computation state struct.
 *
 * Returns: The fraction field of @compstate.
 *
 * Since: 2.49
 **/
gdouble
gwy_computation_state_get_fraction(GwyComputationState *compstate)
{
    return compstate->fraction;
}

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
