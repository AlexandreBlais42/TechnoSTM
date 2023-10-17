/*
 *  $Id: correct.h 21370 2018-08-30 13:34:59Z yeti-dn $
 *  Copyright (C) 2003-2017 David Necas (Yeti), Petr Klapetek.
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

#ifndef __GWY_PROCESS_CORRECT_H__
#define __GWY_PROCESS_CORRECT_H__

#include <libgwyddion/gwymath.h>
#include <libprocess/datafield.h>

G_BEGIN_DECLS

typedef void (*GwyCoordTransform2DFunc)(gdouble x,
                                        gdouble y,
                                        gdouble *px,
                                        gdouble *py,
                                        gpointer user_data);

void         gwy_data_field_laplace_solve               (GwyDataField *field,
                                                         GwyDataField *mask,
                                                         gint grain_id,
                                                         gdouble qprec);
void         gwy_data_field_correct_laplace_iteration   (GwyDataField *data_field,
                                                         GwyDataField *mask_field,
                                                         GwyDataField *buffer_field,
                                                         gdouble corrfactor,
                                                         gdouble *error);
void         gwy_data_field_correct_average             (GwyDataField *data_field,
                                                         GwyDataField *mask_field);
void         gwy_data_field_correct_average_unmasked    (GwyDataField *data_field,
                                                         GwyDataField *mask_field);
void         gwy_data_field_mask_outliers               (GwyDataField *data_field,
                                                         GwyDataField *mask_field,
                                                         gdouble thresh);
void         gwy_data_field_mask_outliers2              (GwyDataField *data_field,
                                                         GwyDataField *mask_field,
                                                         gdouble thresh_low,
                                                         gdouble thresh_high);
void         gwy_data_field_distort                     (GwyDataField *source,
                                                         GwyDataField *dest,
                                                         GwyCoordTransform2DFunc invtrans,
                                                         gpointer user_data,
                                                         GwyInterpolationType interp,
                                                         GwyExteriorType exterior,
                                                         gdouble fill_value);
void         gwy_data_field_sample_distorted            (GwyDataField *source,
                                                         GwyDataField *dest,
                                                         const GwyXY *coords,
                                                         GwyInterpolationType interp,
                                                         GwyExteriorType exterior,
                                                         gdouble fill_value);
void         gwy_data_field_affine                      (GwyDataField *source,
                                                         GwyDataField *dest,
                                                         const gdouble *invtrans,
                                                         GwyInterpolationType interp,
                                                         GwyExteriorType exterior,
                                                         gdouble fill_value);
void         gwy_data_field_affine_prepare              (GwyDataField *source,
                                                         GwyDataField *dest,
                                                         const gdouble *a1a2,
                                                         gdouble *a1a2_corr,
                                                         gdouble *invtrans,
                                                         GwyAffineScalingType scaling,
                                                         gboolean prevent_rotation,
                                                         gdouble oversampling);
gboolean     gwy_data_field_measure_lattice_acf         (GwyDataField *acf2d,
                                                         gdouble *a1a2);
gboolean     gwy_data_field_measure_lattice_psdf        (GwyDataField *psdf2d,
                                                         gdouble *a1a2);
gboolean     gwy_data_line_correct_laplace              (GwyDataLine *data_line,
                                                         GwyDataLine *mask_line);
void         gwy_data_field_mark_scars                  (GwyDataField *data_field,
                                                         GwyDataField *result,
                                                         gdouble threshold_high,
                                                         gdouble threshold_low,
                                                         gdouble min_scar_len,
                                                         gdouble max_scar_width,
                                                         gboolean negative);
void         gwy_data_field_subtract_row_shifts         (GwyDataField *data_field,
                                                         GwyDataLine *shifts);
GwyDataLine* gwy_data_field_find_row_shifts_trimmed_mean(GwyDataField *data_field,
                                                         GwyDataField *mask,
                                                         GwyMaskingType masking,
                                                         gdouble trimfrac,
                                                         gint mincount);
GwyDataLine* gwy_data_field_find_row_shifts_trimmed_diff(GwyDataField *data_field,
                                                         GwyDataField *mask,
                                                         GwyMaskingType masking,
                                                         gdouble trimfrac,
                                                         gint mincount);

GwyPlaneSymmetry gwy_data_field_unrotate_find_corrections(GwyDataLine *derdist,
                                                          gdouble *correction);

G_END_DECLS

#endif /* __GWY_PROCESS_CORRECT__ */

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */
