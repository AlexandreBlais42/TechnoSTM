/*
 *  $Id: filters.h 25188 2023-01-05 16:34:34Z yeti-dn $
 *  Copyright (C) 2003-2019 David Necas (Yeti), Petr Klapetek.
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

#ifndef __GWY_PROCESS_FILTERS_H__
#define __GWY_PROCESS_FILTERS_H__

#include <libprocess/datafield.h>
#include <libprocess/gwyprocessenums.h>

G_BEGIN_DECLS

void     gwy_data_field_normalize                         (GwyDataField *data_field);
void     gwy_data_field_renormalize                       (GwyDataField *data_field,
                                                           gdouble range,
                                                           gdouble offset);
void     gwy_data_field_area_renormalize                  (GwyDataField *data_field,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height,
                                                           gdouble range,
                                                           gdouble offset);
gint     gwy_data_field_threshold                         (GwyDataField *data_field,
                                                           gdouble threshval,
                                                           gdouble bottom,
                                                           gdouble top);
gint     gwy_data_field_area_threshold                    (GwyDataField *data_field,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height,
                                                           gdouble threshval,
                                                           gdouble bottom,
                                                           gdouble top);
gint     gwy_data_field_clamp                             (GwyDataField *data_field,
                                                           gdouble bottom,
                                                           gdouble top);
gint     gwy_data_field_area_clamp                        (GwyDataField *data_field,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height,
                                                           gdouble bottom,
                                                           gdouble top);
void     gwy_data_field_area_gather                       (GwyDataField *data_field,
                                                           GwyDataField *result,
                                                           GwyDataField *buffer,
                                                           gint hsize,
                                                           gint vsize,
                                                           gboolean average,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_convolve                          (GwyDataField *data_field,
                                                           GwyDataField *kernel_field);
void     gwy_data_field_area_convolve                     (GwyDataField *data_field,
                                                           GwyDataField *kernel_field,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_fft_convolve                      (GwyDataField *data_field,
                                                           GwyDataField *kernel_field);
void     gwy_data_field_area_ext_convolve                 (GwyDataField *field,
                                                           guint col,
                                                           guint row,
                                                           guint width,
                                                           guint height,
                                                           GwyDataField *target,
                                                           GwyDataField *kernel,
                                                           GwyExteriorType exterior,
                                                           gdouble fill_value,
                                                           gboolean as_integral);
void     gwy_data_field_convolve_1d                       (GwyDataField *data_field,
                                                           GwyDataLine *kernel_line,
                                                           GwyOrientation orientation);
void     gwy_data_field_area_convolve_1d                  (GwyDataField *data_field,
                                                           GwyDataLine *kernel_line,
                                                           GwyOrientation orientation,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_area_ext_row_convolve             (GwyDataField *field,
                                                           guint col,
                                                           guint row,
                                                           guint width,
                                                           guint height,
                                                           GwyDataField *target,
                                                           GwyDataLine *kernel,
                                                           GwyExteriorType exterior,
                                                           gdouble fill_value,
                                                           gboolean as_integral);
void     gwy_data_field_filter_median                     (GwyDataField *data_field,
                                                           gint size);
void     gwy_data_field_area_filter_median                (GwyDataField *data_field,
                                                           gint size,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_mean                       (GwyDataField *data_field,
                                                           gint size);
void     gwy_data_field_area_filter_mean                  (GwyDataField *data_field,
                                                           gint size,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_conservative               (GwyDataField *data_field,
                                                           gint size);
void     gwy_data_field_area_filter_conservative          (GwyDataField *data_field,
                                                           gint size,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_laplacian                  (GwyDataField *data_field);
void     gwy_data_field_area_filter_laplacian             (GwyDataField *data_field,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_laplacian_of_gaussians     (GwyDataField *data_field);
void     gwy_data_field_area_filter_laplacian_of_gaussians(GwyDataField *data_field,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_sobel                      (GwyDataField *data_field,
                                                           GwyOrientation orientation);
void     gwy_data_field_area_filter_sobel                 (GwyDataField *data_field,
                                                           GwyOrientation orientation,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_sobel_total                (GwyDataField *data_field);
void     gwy_data_field_filter_prewitt                    (GwyDataField *data_field,
                                                           GwyOrientation orientation);
void     gwy_data_field_area_filter_prewitt               (GwyDataField *data_field,
                                                           GwyOrientation orientation,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_prewitt_total              (GwyDataField *data_field);
void     gwy_data_field_filter_slope                      (GwyDataField *data_field,
                                                           GwyDataField *xder,
                                                           GwyDataField *yder);
void     gwy_data_field_filter_gauss_step                 (GwyDataField *data_field,
                                                           gdouble sigma);
void     gwy_data_field_filter_dechecker                  (GwyDataField *data_field);
void     gwy_data_field_area_filter_dechecker             (GwyDataField *data_field,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_gaussian                   (GwyDataField *data_field,
                                                           gdouble sigma);
void     gwy_data_field_area_filter_gaussian              (GwyDataField *data_field,
                                                           gdouble sigma,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_row_gaussian                      (GwyDataField *data_field,
                                                           gdouble sigma);
void     gwy_data_field_column_gaussian                   (GwyDataField *data_field,
                                                           gdouble sigma);
void     gwy_data_field_filter_minimum                    (GwyDataField *data_field,
                                                           gint size);
void     gwy_data_field_area_filter_minimum               (GwyDataField *data_field,
                                                           gint size,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_maximum                    (GwyDataField *data_field,
                                                           gint size);
void     gwy_data_field_area_filter_maximum               (GwyDataField *data_field,
                                                           gint size,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_area_filter_min_max               (GwyDataField *data_field,
                                                           GwyDataField *kernel,
                                                           GwyMinMaxFilterType filtertype,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_area_filter_disc_asf              (GwyDataField *data_field,
                                                           gint radius,
                                                           gboolean closing,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
gboolean gwy_data_field_area_filter_kth_rank              (GwyDataField *data_field,
                                                           GwyDataField *kernel,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height,
                                                           gint k,
                                                           GwySetFractionFunc set_fraction);
void     gwy_data_line_part_filter_kth_rank               (GwyDataLine *data_line,
                                                           gint klen,
                                                           gint from,
                                                           gint to,
                                                           gint k);
gboolean gwy_data_field_area_filter_trimmed_mean          (GwyDataField *data_field,
                                                           GwyDataField *kernel,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height,
                                                           gint nlowest,
                                                           gint nhighest,
                                                           GwySetFractionFunc set_fraction);
void     gwy_data_field_filter_rms                        (GwyDataField *data_field,
                                                           gint size);
void     gwy_data_field_area_filter_rms                   (GwyDataField *data_field,
                                                           gint size,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_kuwahara                   (GwyDataField *data_field);
void     gwy_data_field_area_filter_kuwahara              (GwyDataField *data_field,
                                                           gint col,
                                                           gint row,
                                                           gint width,
                                                           gint height);
void     gwy_data_field_filter_canny                      (GwyDataField *data_field,
                                                           gdouble threshold);
void     gwy_data_field_shade                             (GwyDataField *data_field,
                                                           GwyDataField *target_field,
                                                           gdouble theta,
                                                           gdouble phi);
void     gwy_data_field_filter_harris                     (GwyDataField *x_gradient,
                                                           GwyDataField *y_gradient,
                                                           GwyDataField *result,
                                                           gint neighbourhood,
                                                           gdouble alpha);
void     gwy_data_field_deconvolve_regularized            (GwyDataField *dfield,
                                                           GwyDataField *operand,
                                                           GwyDataField *out,
                                                           gdouble sigma);
void     gwy_data_field_deconvolve_psf_leastsq            (GwyDataField *dfield,
                                                           GwyDataField *operand,
                                                           GwyDataField *out,
                                                           gdouble sigma,
                                                           gint border);
gdouble  gwy_data_field_find_regularization_sigma_for_psf (GwyDataField *dfield,
                                                           GwyDataField *ideal);
gdouble  gwy_data_field_find_regularization_sigma_leastsq (GwyDataField *dfield,
                                                           GwyDataField *ideal,
                                                           gint width,
                                                           gint height,
                                                           gint border);

G_END_DECLS

#endif /* __GWY_PROCESS_FILTERS__ */

/* vim: set cin columns=120 tw=118 et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
