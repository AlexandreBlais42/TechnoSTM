/*
 *  $Id: stats.h 23046 2021-01-08 12:21:05Z yeti-dn $
 *  Copyright (C) 2003-2018 David Necas (Yeti), Petr Klapetek.
 *  E-mail: yeti@gwyddion.net, klapetek@gwyddion.net.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA 02110-1301, USA.
 */

#ifndef __GWY_PROCESS_STATS_H__
#define __GWY_PROCESS_STATS_H__

#include <libprocess/datafield.h>

G_BEGIN_DECLS

gdouble gwy_data_field_get_max              (GwyDataField *data_field);
gdouble gwy_data_field_get_min              (GwyDataField *data_field);
void    gwy_data_field_get_min_max          (GwyDataField *data_field,
                                             gdouble *min,
                                             gdouble *max);
gdouble gwy_data_field_get_avg              (GwyDataField *data_field);
gdouble gwy_data_field_get_rms              (GwyDataField *data_field);
gdouble gwy_data_field_get_mean_square      (GwyDataField *data_field);
gdouble gwy_data_field_get_sum              (GwyDataField *data_field);
gdouble gwy_data_field_get_median           (GwyDataField *data_field);
gdouble gwy_data_field_get_surface_area     (GwyDataField *data_field);
gdouble gwy_data_field_get_surface_slope    (GwyDataField *data_field);
gdouble gwy_data_field_get_variation        (GwyDataField *data_field);
gdouble gwy_data_field_get_entropy          (GwyDataField *data_field);
gdouble gwy_data_field_get_entropy_2d       (GwyDataField *xfield,
                                             GwyDataField *yfield);
gdouble gwy_data_field_area_get_max         (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
gdouble gwy_data_field_area_get_min         (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
#ifndef GWY_DISABLE_DEPRECATED
void    gwy_data_field_area_get_min_max     (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height,
                                             gdouble *min,
                                             gdouble *max);
#endif
void    gwy_data_field_area_get_min_max_mask(GwyDataField *data_field,
                                             GwyDataField *mask,
                                             GwyMaskingType mode,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height,
                                             gdouble *min,
                                             gdouble *max);
#ifndef GWY_DISABLE_DEPRECATED
gdouble gwy_data_field_area_get_avg         (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
#endif
gdouble gwy_data_field_area_get_avg_mask    (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             GwyMaskingType mode,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
#ifndef GWY_DISABLE_DEPRECATED
gdouble gwy_data_field_area_get_rms         (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
#endif
gdouble gwy_data_field_area_get_rms_mask    (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             GwyMaskingType mode,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
gdouble gwy_data_field_area_get_grainwise_rms(GwyDataField *data_field,
                                              GwyDataField *mask,
                                              GwyMaskingType mode,
                                              gint col,
                                              gint row,
                                              gint width,
                                              gint height);
#ifndef GWY_DISABLE_DEPRECATED
gdouble gwy_data_field_area_get_sum         (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
#endif
gdouble gwy_data_field_area_get_sum_mask    (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             GwyMaskingType mode,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
#ifndef GWY_DISABLE_DEPRECATED
gdouble gwy_data_field_area_get_median      (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
#endif
gdouble gwy_data_field_area_get_median_mask (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             GwyMaskingType mode,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
#ifndef GWY_DISABLE_DEPRECATED
gdouble gwy_data_field_area_get_surface_area(GwyDataField *data_field,
                                             GwyDataField *mask,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
#endif
gdouble gwy_data_field_area_get_surface_area_mask(GwyDataField *data_field,
                                                  GwyDataField *mask,
                                                  GwyMaskingType mode,
                                                  gint col,
                                                  gint row,
                                                  gint width,
                                                  gint height);
gdouble gwy_data_field_area_get_surface_slope_mask(GwyDataField *data_field,
                                                   GwyDataField *mask,
                                                   GwyMaskingType mode,
                                                   gint col,
                                                   gint row,
                                                   gint width,
                                                   gint height);
gdouble gwy_data_field_area_get_mean_square (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             GwyMaskingType mode,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
gdouble gwy_data_field_area_get_entropy_at_scales(GwyDataField *data_field,
                                                  GwyDataLine *target_line,
                                                  GwyDataField *mask,
                                                  GwyMaskingType mode,
                                                  gint col,
                                                  gint row,
                                                  gint width,
                                                  gint height,
                                                  gint maxdiv);
gdouble gwy_data_field_get_entropy_2d_at_scales  (GwyDataField *xfield,
                                                  GwyDataField *yfield,
                                                  GwyDataLine *target_line,
                                                  gint maxdiv);
gdouble gwy_data_field_area_get_variation   (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             GwyMaskingType mode,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
gdouble gwy_data_field_area_get_entropy     (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             GwyMaskingType mode,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
gdouble gwy_data_field_area_get_volume      (GwyDataField *data_field,
                                             GwyDataField *basis,
                                             GwyDataField *mask,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height);
void    gwy_data_field_get_autorange        (GwyDataField *data_field,
                                             gdouble *from,
                                             gdouble *to);
void    gwy_data_field_area_get_autorange   (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             GwyMaskingType masking,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height,
                                             gdouble *from,
                                             gdouble *to);
void    gwy_data_field_get_stats            (GwyDataField *data_field,
                                             gdouble *avg,
                                             gdouble *ra,
                                             gdouble *rms,
                                             gdouble *skew,
                                             gdouble *kurtosis);

#ifndef GWY_DISABLE_DEPRECATED
void    gwy_data_field_area_get_stats       (GwyDataField *data_field,
                                             GwyDataField *mask,
                                             gint col,
                                             gint row,
                                             gint width,
                                             gint height,
                                             gdouble *avg,
                                             gdouble *ra,
                                             gdouble *rms,
                                             gdouble *skew,
                                             gdouble *kurtosis);
#endif

void         gwy_data_field_area_get_stats_mask     (GwyDataField *data_field,
                                                     GwyDataField *mask,
                                                     GwyMaskingType mode,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gdouble *avg,
                                                     gdouble *ra,
                                                     gdouble *rms,
                                                     gdouble *skew,
                                                     gdouble *kurtosis);
void         gwy_data_field_area_count_in_range     (GwyDataField *data_field,
                                                     GwyDataField *mask,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gdouble below,
                                                     gdouble above,
                                                     gint *nbelow,
                                                     gint *nabove);
void         gwy_data_field_area_dh                 (GwyDataField *data_field,
                                                     GwyDataField *mask,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gint nstats);
void         gwy_data_field_dh                      (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint nstats);
void         gwy_data_field_area_cdh                (GwyDataField *data_field,
                                                     GwyDataField *mask,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gint nstats);
void         gwy_data_field_cdh                     (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint nstats);
void         gwy_data_field_area_da                 (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyOrientation orientation,
                                                     gint nstats);
void         gwy_data_field_area_da_mask            (GwyDataField *data_field,
                                                     GwyDataField *mask,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyOrientation orientation,
                                                     gint nstats);
void         gwy_data_field_da                      (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     GwyOrientation orientation,
                                                     gint nstats);
void         gwy_data_field_area_cda                (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyOrientation orientation,
                                                     gint nstats);
void         gwy_data_field_area_cda_mask           (GwyDataField *data_field,
                                                     GwyDataField *mask,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyOrientation orientation,
                                                     gint nstats);
void         gwy_data_field_cda                     (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     GwyOrientation orientation,
                                                     gint nstats);
void         gwy_data_field_area_acf                (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyOrientation orientation,
                                                     GwyInterpolationType interpolation,
                                                     gint nstats);
void         gwy_data_field_acf                     (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     GwyOrientation orientation,
                                                     GwyInterpolationType interpolation,
                                                     gint nstats);
GwyDataLine* gwy_data_field_area_row_acf            (GwyDataField *field,
                                                     GwyDataField *mask,
                                                     GwyMaskingType masking,
                                                     guint col,
                                                     guint row,
                                                     guint width,
                                                     guint height,
                                                     guint level,
                                                     GwyDataLine *weights);
void         gwy_data_field_area_hhcf               (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyOrientation orientation,
                                                     GwyInterpolationType interpolation,
                                                     gint nstats);
void         gwy_data_field_hhcf                    (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     GwyOrientation orientation,
                                                     GwyInterpolationType interpolation,
                                                     gint nstats);
GwyDataLine* gwy_data_field_area_row_hhcf           (GwyDataField *field,
                                                     GwyDataField *mask,
                                                     GwyMaskingType masking,
                                                     guint col,
                                                     guint row,
                                                     guint width,
                                                     guint height,
                                                     guint level,
                                                     GwyDataLine *weights);
void         gwy_data_field_area_psdf               (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyOrientation orientation,
                                                     GwyInterpolationType interpolation,
                                                     GwyWindowingType windowing,
                                                     gint nstats);
void         gwy_data_field_psdf                    (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     GwyOrientation orientation,
                                                     GwyInterpolationType interpolation,
                                                     GwyWindowingType windowing,
                                                     gint nstats);
GwyDataLine* gwy_data_field_area_row_psdf           (GwyDataField *field,
                                                     GwyDataField *mask,
                                                     GwyMaskingType masking,
                                                     guint col,
                                                     guint row,
                                                     guint width,
                                                     guint height,
                                                     GwyWindowingType windowing,
                                                     guint level);
void         gwy_data_field_area_rpsdf              (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyInterpolationType interpolation,
                                                     GwyWindowingType windowing,
                                                     gint nstats);
void         gwy_data_field_rpsdf                   (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     GwyInterpolationType interpolation,
                                                     GwyWindowingType windowing,
                                                     gint nstats);
GwyDataLine* gwy_data_field_area_row_asg            (GwyDataField *field,
                                                     GwyDataField *mask,
                                                     GwyMaskingType masking,
                                                     guint col,
                                                     guint row,
                                                     guint width,
                                                     guint height,
                                                     guint level);
void         gwy_data_field_area_2dacf              (GwyDataField *data_field,
                                                     GwyDataField *target_field,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gint xrange,
                                                     gint yrange);
void         gwy_data_field_area_2dacf_mask         (GwyDataField *data_field,
                                                     GwyDataField *target_field,
                                                     GwyDataField *mask,
                                                     GwyMaskingType masking,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gint xrange,
                                                     gint yrange,
                                                     GwyDataField *weights);
void         gwy_data_field_2dacf                   (GwyDataField *data_field,
                                                     GwyDataField *target_field);
void         gwy_data_field_area_2dpsdf_mask        (GwyDataField *field,
                                                     GwyDataField *target_field,
                                                     GwyDataField *mask,
                                                     GwyMaskingType masking,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyWindowingType windowing,
                                                     gint level);
void         gwy_data_field_2dpsdf                  (GwyDataField *field,
                                                     GwyDataField *target_field,
                                                     GwyWindowingType windowing,
                                                     gint level);
void         gwy_data_field_area_racf               (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gint nstats);
void         gwy_data_field_racf                    (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint nstats);
void         gwy_data_field_area_minkowski_volume   (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gint nstats);
void         gwy_data_field_minkowski_volume        (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint nstats);
void         gwy_data_field_area_minkowski_boundary (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gint nstats);
void         gwy_data_field_minkowski_boundary      (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint nstats);
void         gwy_data_field_area_minkowski_euler    (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gint nstats);
void         gwy_data_field_minkowski_euler         (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     gint nstats);
gdouble      gwy_data_field_area_get_dispersion     (GwyDataField *data_field,
                                                     GwyDataField *mask,
                                                     GwyMaskingType masking,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gdouble *xcenter,
                                                     gdouble *ycenter);
gdouble      gwy_data_field_get_dispersion          (GwyDataField *data_field,
                                                     gdouble *xcenter,
                                                     gdouble *ycenter);
void         gwy_data_field_slope_distribution      (GwyDataField *data_field,
                                                     GwyDataLine *derdist,
                                                     gint kernel_size);
void         gwy_data_field_get_normal_coeffs       (GwyDataField *data_field,
                                                     gdouble *nx,
                                                     gdouble *ny,
                                                     gdouble *nz,
                                                     gboolean normalize1);
void         gwy_data_field_area_get_normal_coeffs  (GwyDataField *data_field,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gdouble *nx,
                                                     gdouble *ny,
                                                     gdouble *nz,
                                                     gboolean normalize1);
void         gwy_data_field_area_get_inclination    (GwyDataField *data_field,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     gdouble *theta,
                                                     gdouble *phi);
void         gwy_data_field_get_inclination         (GwyDataField *data_field,
                                                     gdouble *theta,
                                                     gdouble *phi);
void         gwy_data_field_area_get_line_stats     (GwyDataField *data_field,
                                                     GwyDataField *mask,
                                                     GwyDataLine *target_line,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyLineStatQuantity quantity,
                                                     GwyOrientation orientation);
void         gwy_data_field_get_line_stats_mask     (GwyDataField *data_field,
                                                     GwyDataField *mask,
                                                     GwyMaskingType masking,
                                                     GwyDataLine *target_line,
                                                     GwyDataLine *weights,
                                                     gint col,
                                                     gint row,
                                                     gint width,
                                                     gint height,
                                                     GwyLineStatQuantity quantity,
                                                     GwyOrientation orientation);
void         gwy_data_field_get_line_stats          (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     GwyLineStatQuantity quantity,
                                                     GwyOrientation orientation);
guint        gwy_data_field_count_maxima            (GwyDataField *data_field);
guint        gwy_data_field_count_minima            (GwyDataField *data_field);
GwyDataLine* gwy_data_field_psdf_to_angular_spectrum(GwyDataField *psdf,
                                                     gint nstats);
void         gwy_data_field_angular_average         (GwyDataField *data_field,
                                                     GwyDataLine *target_line,
                                                     GwyDataField *mask,
                                                     GwyMaskingType masking,
                                                     gdouble x,
                                                     gdouble y,
                                                     gdouble r,
                                                     gint nstats);

G_END_DECLS

#endif /* __GWY_PROCESS_STATS_H__ */

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

