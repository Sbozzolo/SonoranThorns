! Copyright (C) 2021 Gabriele Bozzola

! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 3 of the License, or (at your option) any later
! version.

! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.

! You should have received a copy of the GNU General Public License along with
! this program; if not, see <https://www.gnu.org/licenses/>.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Christoffel_zero( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  Gamma_ttt = 0
  Gamma_ttx = 0
  Gamma_tty = 0
  Gamma_ttz = 0
  Gamma_txx = 0
  Gamma_txy = 0
  Gamma_txz = 0
  Gamma_tyy = 0
  Gamma_tyz = 0
  Gamma_tzz = 0
  Gamma_xtt = 0
  Gamma_xtx = 0
  Gamma_xty = 0
  Gamma_xtz = 0
  Gamma_xxx = 0
  Gamma_xxy = 0
  Gamma_xxz = 0
  Gamma_xyy = 0
  Gamma_xyz = 0
  Gamma_xzz = 0
  Gamma_ytt = 0
  Gamma_ytx = 0
  Gamma_yty = 0
  Gamma_ytz = 0
  Gamma_yxx = 0
  Gamma_yxy = 0
  Gamma_yxz = 0
  Gamma_yyy = 0
  Gamma_yyz = 0
  Gamma_yzz = 0
  Gamma_ztt = 0
  Gamma_ztx = 0
  Gamma_zty = 0
  Gamma_ztz = 0
  Gamma_zxx = 0
  Gamma_zxy = 0
  Gamma_zxz = 0
  Gamma_zyy = 0
  Gamma_zyz = 0
  Gamma_zzz = 0

  if (save_dgab/=0) then
     g_ttt = 0
     g_ttx = 0
     g_tty = 0
     g_ttz = 0
     g_txx = 0
     g_txy = 0
     g_txz = 0
     g_tyy = 0
     g_tyz = 0
     g_tzz = 0
     g_xtt = 0
     g_xtx = 0
     g_xty = 0
     g_xtz = 0
     g_xxx = 0
     g_xxy = 0
     g_xxz = 0
     g_xyy = 0
     g_xyz = 0
     g_xzz = 0
     g_ytt = 0
     g_ytx = 0
     g_yty = 0
     g_ytz = 0
     g_yxx = 0
     g_yxy = 0
     g_yxz = 0
     g_yyy = 0
     g_yyz = 0
     g_yzz = 0
     g_ztt = 0
     g_ztx = 0
     g_zty = 0
     g_ztz = 0
     g_zxx = 0
     g_zxy = 0
     g_zxz = 0
     g_zyy = 0
     g_zyz = 0
     g_zzz = 0
  end if

end subroutine Christoffel_zero

subroutine Christoffel_compute( CCTK_ARGUMENTS )

  use adm_metric

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT i, j, k
  CCTK_INT a, b, c, m
  CCTK_INT ierr

  ! Here we are not using the symmetry of the matrix
  CCTK_REAL :: gg(3,3), ggu(3,3), gab(4,4), gUP(4,4), dgg(3,3,3), dgab(4,4,4)
  CCTK_REAL :: cf1(4,4,4), cf2(4,4,4)
  CCTK_REAL :: alph, one_over_alph_squared, detg, beta(3)
  CCTK_REAL :: dt_alph, dt_beta(3), dt_gg(3,3)
  CCTK_REAL :: dalph(3), dbeta(3,3)
  CCTK_REAL :: ww, dww, oww_sq
  CCTK_REAL, parameter ::  one  = 1.d0

  CCTK_REAL dx2, dy2, dz2, dx12, dy12, dz12, odx60, ody60, odz60

  if (compute_every .le. 0) then
     goto 9999
  end if

  if (MOD(cctk_iteration, compute_every) .ne. 0) then
     goto 9999
  endif

  dx2 = 2*CCTK_DELTA_SPACE(1)
  dy2 = 2*CCTK_DELTA_SPACE(2)
  dz2 = 2*CCTK_DELTA_SPACE(3)

  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)

  odx60 = 1 / (60 * CCTK_DELTA_SPACE(1))
  ody60 = 1 / (60 * CCTK_DELTA_SPACE(2))
  odz60 = 1 / (60 * CCTK_DELTA_SPACE(3))

  !!$OMP PARALLEL DO COLLAPSE(3) &
  !!$OMP PRIVATE( i, j, k, gab, dgab, ww, dww, oww_sq, &
  !!$OMP gg, ggu, alph, one_over_alph_squared, detg, &
  !!$OMP beta, dt_alph, dt_beta, dt_gg,  &
  !!$OMP dgg, dalph, dbeta, cf1, cf2)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
     do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
        do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

           gg(1,1) = gxx(i,j,k)
           gg(1,2) = gxy(i,j,k)
           gg(1,3) = gxz(i,j,k)
           gg(2,2) = gyy(i,j,k)
           gg(2,3) = gyz(i,j,k)
           gg(3,3) = gzz(i,j,k)
           gg(2,1) = gg(1,2)
           gg(3,1) = gg(1,3)
           gg(3,2) = gg(2,3)

           alph      = alp(i,j,k)

           one_over_alph_squared = one / (alph * alph)

           beta(1)   = betax(i,j,k)
           beta(2)   = betay(i,j,k)
           beta(3)   = betaz(i,j,k)

           dt_alph      = dtalp(i,j,k)
           dt_beta(1)   = dtbetax(i,j,k)
           dt_beta(2)   = dtbetay(i,j,k)
           dt_beta(3)   = dtbetaz(i,j,k)

           ww = conf_fac(i,j,k)

           ! conf_fac = detg**-1/6
           detg = ww**(-6.d0)

           dww = - 2.d0 * rhs_conf_fac(i, j, k) / (ww * ww * ww)
           oww_sq = one / (ww * ww)

           dt_gg(1,1) = oww_sq * rhs_hxx(i,j,k) + dww * hxx(i,j,k)
           dt_gg(1,2) = oww_sq * rhs_hxy(i,j,k) + dww * hxy(i,j,k)
           dt_gg(1,3) = oww_sq * rhs_hxz(i,j,k) + dww * hxz(i,j,k)
           dt_gg(2,2) = oww_sq * rhs_hyy(i,j,k) + dww * hyy(i,j,k)
           dt_gg(2,3) = oww_sq * rhs_hyz(i,j,k) + dww * hyz(i,j,k)
           dt_gg(3,3) = oww_sq * rhs_hzz(i,j,k) + dww * hzz(i,j,k)
           dt_gg(2,1) = dt_gg(1,2)
           dt_gg(3,1) = dt_gg(1,3)
           dt_gg(3,2) = dt_gg(2,3)

           call calc_4metric(gg, alph, beta, gab)

           if (derivs_order == 6) then

              dgg(1,1,1) = (  gxx(i+3,j,k) - 9*gxx(i+2,j,k) + 45*gxx(i+1,j,k)      &
                   - gxx(i-3,j,k) + 9*gxx(i-2,j,k) - 45*gxx(i-1,j,k) ) * odx60
              dgg(1,2,1) = (  gxy(i+3,j,k) - 9*gxy(i+2,j,k) + 45*gxy(i+1,j,k)      &
                   - gxy(i-3,j,k) + 9*gxy(i-2,j,k) - 45*gxy(i-1,j,k) ) * odx60
              dgg(1,3,1) = (  gxz(i+3,j,k) - 9*gxz(i+2,j,k) + 45*gxz(i+1,j,k)      &
                   - gxz(i-3,j,k) + 9*gxz(i-2,j,k) - 45*gxz(i-1,j,k) ) * odx60
              dgg(2,2,1) = (  gyy(i+3,j,k) - 9*gyy(i+2,j,k) + 45*gyy(i+1,j,k)      &
                   - gyy(i-3,j,k) + 9*gyy(i-2,j,k) - 45*gyy(i-1,j,k) ) * odx60
              dgg(2,3,1) = (  gyz(i+3,j,k) - 9*gyz(i+2,j,k) + 45*gyz(i+1,j,k)      &
                   - gyz(i-3,j,k) + 9*gyz(i-2,j,k) - 45*gyz(i-1,j,k) ) * odx60
              dgg(3,3,1) = (  gzz(i+3,j,k) - 9*gzz(i+2,j,k) + 45*gzz(i+1,j,k)      &
                   - gzz(i-3,j,k) + 9*gzz(i-2,j,k) - 45*gzz(i-1,j,k) ) * odx60

              dgg(1,1,2) = (  gxx(i,j+3,k) - 9*gxx(i,j+2,k) + 45*gxx(i,j+1,k)      &
                   - gxx(i,j-3,k) + 9*gxx(i,j-2,k) - 45*gxx(i,j-1,k) ) * ody60
              dgg(1,2,2) = (  gxy(i,j+3,k) - 9*gxy(i,j+2,k) + 45*gxy(i,j+1,k)      &
                   - gxy(i,j-3,k) + 9*gxy(i,j-2,k) - 45*gxy(i,j-1,k) ) * ody60
              dgg(1,3,2) = (  gxz(i,j+3,k) - 9*gxz(i,j+2,k) + 45*gxz(i,j+1,k)      &
                   - gxz(i,j-3,k) + 9*gxz(i,j-2,k) - 45*gxz(i,j-1,k) ) * ody60
              dgg(2,2,2) = (  gyy(i,j+3,k) - 9*gyy(i,j+2,k) + 45*gyy(i,j+1,k)      &
                   - gyy(i,j-3,k) + 9*gyy(i,j-2,k) - 45*gyy(i,j-1,k) ) * ody60
              dgg(2,3,2) = (  gyz(i,j+3,k) - 9*gyz(i,j+2,k) + 45*gyz(i,j+1,k)      &
                   - gyz(i,j-3,k) + 9*gyz(i,j-2,k) - 45*gyz(i,j-1,k) ) * ody60
              dgg(3,3,2) = (  gzz(i,j+3,k) - 9*gzz(i,j+2,k) + 45*gzz(i,j+1,k)      &
                   - gzz(i,j-3,k) + 9*gzz(i,j-2,k) - 45*gzz(i,j-1,k) ) * ody60

              dgg(1,1,3) = (  gxx(i,j,k+3) - 9*gxx(i,j,k+2) + 45*gxx(i,j,k+1)      &
                   - gxx(i,j,k-3) + 9*gxx(i,j,k-2) - 45*gxx(i,j,k-1) ) * odz60
              dgg(1,2,3) = (  gxy(i,j,k+3) - 9*gxy(i,j,k+2) + 45*gxy(i,j,k+1)      &
                   - gxy(i,j,k-3) + 9*gxy(i,j,k-2) - 45*gxy(i,j,k-1) ) * odz60
              dgg(1,3,3) = (  gxz(i,j,k+3) - 9*gxz(i,j,k+2) + 45*gxz(i,j,k+1)      &
                   - gxz(i,j,k-3) + 9*gxz(i,j,k-2) - 45*gxz(i,j,k-1) ) * odz60
              dgg(2,2,3) = (  gyy(i,j,k+3) - 9*gyy(i,j,k+2) + 45*gyy(i,j,k+1)      &
                   - gyy(i,j,k-3) + 9*gyy(i,j,k-2) - 45*gyy(i,j,k-1) ) * odz60
              dgg(2,3,3) = (  gyz(i,j,k+3) - 9*gyz(i,j,k+2) + 45*gyz(i,j,k+1)      &
                   - gyz(i,j,k-3) + 9*gyz(i,j,k-2) - 45*gyz(i,j,k-1) ) * odz60
              dgg(3,3,3) = (  gzz(i,j,k+3) - 9*gzz(i,j,k+2) + 45*gzz(i,j,k+1)      &
                   - gzz(i,j,k-3) + 9*gzz(i,j,k-2) - 45*gzz(i,j,k-1) ) * odz60

              dgg(2,1,:) = dgg(1,2,:)
              dgg(3,1,:) = dgg(1,3,:)
              dgg(3,2,:) = dgg(2,3,:)

              dalph(1) = (  alp(i+3,j,k) - 9*alp(i+2,j,k) + 45*alp(i+1,j,k) &
                   - alp(i-3,j,k) + 9*alp(i-2,j,k) - 45*alp(i-1,j,k) ) * odx60

              dalph(2) = (  alp(i,j+3,k) - 9*alp(i,j+2,k) + 45*alp(i,j+1,k) &
                   - alp(i,j-3,k) + 9*alp(i,j-2,k) - 45*alp(i,j-1,k) ) * ody60

              dalph(3) = (  alp(i,j,k+3) - 9*alp(i,j,k+2) + 45*alp(i,j,k+1) &
                   - alp(i,j,k-3) + 9*alp(i,j,k-2) - 45*alp(i,j,k-1) ) * odz60


              dbeta(1,1) = (  betax(i+3,j,k) - 9*betax(i+2,j,k) + 45*betax(i+1,j,k) &
                   - betax(i-3,j,k) + 9*betax(i-2,j,k) - 45*betax(i-1,j,k) ) * odx60
              dbeta(2,1) = (  betay(i+3,j,k) - 9*betay(i+2,j,k) + 45*betay(i+1,j,k) &
                   - betay(i-3,j,k) + 9*betay(i-2,j,k) - 45*betay(i-1,j,k) ) * odx60
              dbeta(3,1) = (  betaz(i+3,j,k) - 9*betaz(i+2,j,k) + 45*betaz(i+1,j,k) &
                   - betaz(i-3,j,k) + 9*betaz(i-2,j,k) - 45*betaz(i-1,j,k) ) * odx60

              dbeta(1,2) = (  betax(i,j+3,k) - 9*betax(i,j+2,k) + 45*betax(i,j+1,k) &
                   - betax(i,j-3,k) + 9*betax(i,j-2,k) - 45*betax(i,j-1,k) ) * ody60
              dbeta(2,2) = (  betay(i,j+3,k) - 9*betay(i,j+2,k) + 45*betay(i,j+1,k) &
                   - betay(i,j-3,k) + 9*betay(i,j-2,k) - 45*betay(i,j-1,k) ) * ody60
              dbeta(3,2) = (  betaz(i,j+3,k) - 9*betaz(i,j+2,k) + 45*betaz(i,j+1,k) &
                   - betaz(i,j-3,k) + 9*betaz(i,j-2,k) - 45*betaz(i,j-1,k) ) * ody60

              dbeta(1,3) = (  betax(i,j,k+3) - 9*betax(i,j,k+2) + 45*betax(i,j,k+1) &
                   - betax(i,j,k-3) + 9*betax(i,j,k-2) - 45*betax(i,j,k-1) ) * odz60
              dbeta(2,3) = (  betay(i,j,k+3) - 9*betay(i,j,k+2) + 45*betay(i,j,k+1) &
                   - betay(i,j,k-3) + 9*betay(i,j,k-2) - 45*betay(i,j,k-1) ) * odz60
              dbeta(3,3) = (  betaz(i,j,k+3) - 9*betaz(i,j,k+2) + 45*betaz(i,j,k+1) &
                   - betaz(i,j,k-3) + 9*betaz(i,j,k-2) - 45*betaz(i,j,k-1) ) * odz60

           else if (derivs_order == 4) then
              dgg(1,1,1) = (   -gxx(i+2,j,k) + 8*gxx(i+1,j,k)                      &
                   - 8*gxx(i-1,j,k) +   gxx(i-2,j,k) ) / dx12
              dgg(1,2,1) = (   -gxy(i+2,j,k) + 8*gxy(i+1,j,k)                      &
                   - 8*gxy(i-1,j,k) +   gxy(i-2,j,k) ) / dx12
              dgg(1,3,1) = (   -gxz(i+2,j,k) + 8*gxz(i+1,j,k)                      &
                   - 8*gxz(i-1,j,k) +   gxz(i-2,j,k) ) / dx12
              dgg(2,2,1) = (   -gyy(i+2,j,k) + 8*gyy(i+1,j,k)                      &
                   - 8*gyy(i-1,j,k) +   gyy(i-2,j,k) ) / dx12
              dgg(2,3,1) = (   -gyz(i+2,j,k) + 8*gyz(i+1,j,k)                      &
                   - 8*gyz(i-1,j,k) +   gyz(i-2,j,k) ) / dx12
              dgg(3,3,1) = (   -gzz(i+2,j,k) + 8*gzz(i+1,j,k)                      &
                   - 8*gzz(i-1,j,k) +   gzz(i-2,j,k) ) / dx12

              dgg(1,1,2) = (   -gxx(i,j+2,k) + 8*gxx(i,j+1,k)                      &
                   - 8*gxx(i,j-1,k) +   gxx(i,j-2,k) ) / dy12
              dgg(1,2,2) = (   -gxy(i,j+2,k) + 8*gxy(i,j+1,k)                      &
                   - 8*gxy(i,j-1,k) +   gxy(i,j-2,k) ) / dy12
              dgg(1,3,2) = (   -gxz(i,j+2,k) + 8*gxz(i,j+1,k)                      &
                   - 8*gxz(i,j-1,k) +   gxz(i,j-2,k) ) / dy12
              dgg(2,2,2) = (   -gyy(i,j+2,k) + 8*gyy(i,j+1,k)                      &
                   - 8*gyy(i,j-1,k) +   gyy(i,j-2,k) ) / dy12
              dgg(2,3,2) = (   -gyz(i,j+2,k) + 8*gyz(i,j+1,k)                      &
                   - 8*gyz(i,j-1,k) +   gyz(i,j-2,k) ) / dy12
              dgg(3,3,2) = (   -gzz(i,j+2,k) + 8*gzz(i,j+1,k)                      &
                   - 8*gzz(i,j-1,k) +   gzz(i,j-2,k) ) / dy12

              dgg(1,1,3) = (   -gxx(i,j,k+2) + 8*gxx(i,j,k+1)                      &
                   - 8*gxx(i,j,k-1) +   gxx(i,j,k-2) ) / dz12
              dgg(1,2,3) = (   -gxy(i,j,k+2) + 8*gxy(i,j,k+1)                      &
                   - 8*gxy(i,j,k-1) +   gxy(i,j,k-2) ) / dz12
              dgg(1,3,3) = (   -gxz(i,j,k+2) + 8*gxz(i,j,k+1)                      &
                   - 8*gxz(i,j,k-1) +   gxz(i,j,k-2) ) / dz12
              dgg(2,2,3) = (   -gyy(i,j,k+2) + 8*gyy(i,j,k+1)                      &
                   - 8*gyy(i,j,k-1) +   gyy(i,j,k-2) ) / dz12
              dgg(2,3,3) = (   -gyz(i,j,k+2) + 8*gyz(i,j,k+1)                      &
                   - 8*gyz(i,j,k-1) +   gyz(i,j,k-2) ) / dz12
              dgg(3,3,3) = (   -gzz(i,j,k+2) + 8*gzz(i,j,k+1)                      &
                   - 8*gzz(i,j,k-1) +   gzz(i,j,k-2) ) / dz12


              dgg(2,1,:) = dgg(1,2,:)
              dgg(3,1,:) = dgg(1,3,:)
              dgg(3,2,:) = dgg(2,3,:)

              dalph(1) = (   -alp(i+2,j,k) + 8*alp(i+1,j,k)                        &
                   - 8*alp(i-1,j,k) +   alp(i-2,j,k) ) / dx12

              dalph(2) = (   -alp(i,j+2,k) + 8*alp(i,j+1,k)                        &
                   - 8*alp(i,j-1,k) +   alp(i,j-2,k) ) / dy12

              dalph(3) = (   -alp(i,j,k+2) + 8*alp(i,j,k+1)                        &
                   - 8*alp(i,j,k-1) +   alp(i,j,k-2) ) / dz12

              dbeta(1,1)  = (   -betax(i+2,j,k) + 8*betax(i+1,j,k)                 &
                   - 8*betax(i-1,j,k) +   betax(i-2,j,k) ) / dx12
              dbeta(2,1)  = (   -betay(i+2,j,k) + 8*betay(i+1,j,k)                 &
                   - 8*betay(i-1,j,k) +   betay(i-2,j,k) ) / dx12
              dbeta(3,1)  = (   -betaz(i+2,j,k) + 8*betaz(i+1,j,k)                 &
                   - 8*betaz(i-1,j,k) +   betaz(i-2,j,k) ) / dx12

              dbeta(1,2)  = (   -betax(i,j+2,k) + 8*betax(i,j+1,k)                 &
                   - 8*betax(i,j-1,k) +   betax(i,j-2,k) ) / dy12
              dbeta(2,2)  = (   -betay(i,j+2,k) + 8*betay(i,j+1,k)                 &
                   - 8*betay(i,j-1,k) +   betay(i,j-2,k) ) / dy12
              dbeta(3,2)  = (   -betaz(i,j+2,k) + 8*betaz(i,j+1,k)                 &
                   - 8*betaz(i,j-1,k) +   betaz(i,j-2,k) ) / dy12

              dbeta(1,3)  = (   -betax(i,j,k+2) + 8*betax(i,j,k+1)                 &
                   - 8*betax(i,j,k-1) +   betax(i,j,k-2) ) / dz12
              dbeta(2,3)  = (   -betay(i,j,k+2) + 8*betay(i,j,k+1)                 &
                   - 8*betay(i,j,k-1) +   betay(i,j,k-2) ) / dz12
              dbeta(3,3)  = (   -betaz(i,j,k+2) + 8*betaz(i,j,k+1)                 &
                   - 8*betaz(i,j,k-1) +   betaz(i,j,k-2) ) / dz12

           else if (derivs_order == 2) then
              dgg(1,1,1) = (gxx(i+1,j,k) - gxx(i-1,j,k)) / dx2
              dgg(1,2,1) = (gxy(i+1,j,k) - gxy(i-1,j,k)) / dx2
              dgg(1,3,1) = (gxz(i+1,j,k) - gxz(i-1,j,k)) / dx2
              dgg(2,2,1) = (gyy(i+1,j,k) - gyy(i-1,j,k)) / dx2
              dgg(2,3,1) = (gyz(i+1,j,k) - gyz(i-1,j,k)) / dx2
              dgg(3,3,1) = (gzz(i+1,j,k) - gzz(i-1,j,k)) / dx2

              dgg(1,1,2) = (gxx(i,j+1,k) - gxx(i,j-1,k)) / dy2
              dgg(1,2,2) = (gxy(i,j+1,k) - gxy(i,j-1,k)) / dy2
              dgg(1,3,2) = (gxz(i,j+1,k) - gxz(i,j-1,k)) / dy2
              dgg(2,2,2) = (gyy(i,j+1,k) - gyy(i,j-1,k)) / dy2
              dgg(2,3,2) = (gyz(i,j+1,k) - gyz(i,j-1,k)) / dy2
              dgg(3,3,2) = (gzz(i,j+1,k) - gzz(i,j-1,k)) / dy2

              dgg(1,1,3) = (gxx(i,j,k+1) - gxx(i,j,k-1)) / dz2
              dgg(1,2,3) = (gxy(i,j,k+1) - gxy(i,j,k-1)) / dz2
              dgg(1,3,3) = (gxz(i,j,k+1) - gxz(i,j,k-1)) / dz2
              dgg(2,2,3) = (gyy(i,j,k+1) - gyy(i,j,k-1)) / dz2
              dgg(2,3,3) = (gyz(i,j,k+1) - gyz(i,j,k-1)) / dz2
              dgg(3,3,3) = (gzz(i,j,k+1) - gzz(i,j,k-1)) / dz2


              dgg(2,1,:) = dgg(1,2,:)
              dgg(3,1,:) = dgg(1,3,:)
              dgg(3,2,:) = dgg(2,3,:)

              dalph(1) = (alp(i+1,j,k) - alp(i-1,j,k)) / dx2

              dalph(2) = (alp(i,j+1,k) - alp(i,j-1,k)) / dy2

              dalph(3) = (alp(i,j,k+1) - alp(i,j,k-1)) / dz2

              dbeta(1,1)  = (betax(i+1,j,k) - betax(i-1,j,k)) / dx2
              dbeta(2,1)  = (betay(i+1,j,k) - betay(i-1,j,k)) / dx2
              dbeta(3,1)  = (betaz(i+1,j,k) - betaz(i-1,j,k)) / dx2

              dbeta(1,2)  = (betax(i,j+1,k) - betax(i,j-1,k)) / dy2
              dbeta(2,2)  = (betay(i,j+1,k) - betay(i,j-1,k)) / dy2
              dbeta(3,2)  = (betaz(i,j+1,k) - betaz(i,j-1,k)) / dy2

              dbeta(1,3)  = (betax(i,j,k+1) - betax(i,j,k-1)) / dz2
              dbeta(2,3)  = (betay(i,j,k+1) - betay(i,j,k-1)) / dz2
              dbeta(3,3)  = (betaz(i,j,k+1) - betaz(i,j,k-1)) / dz2
           else
              call CCTK_ERROR("derivs_order not implemented")
           end if

           call calc_4metricderivs(gg, alph, beta, dgg, dalph, dbeta, &
                dt_gg, dt_alph, dt_beta, gab, dgab)

           ! Inverse
           detg    =       gg(1,1) * gg(2,2) * gg(3,3)                            &
                + 2 * gg(1,2) * gg(1,3) * gg(2,3)                            &
                -     gg(1,1) * gg(2,3) ** 2                                 &
                -     gg(2,2) * gg(1,3) ** 2                                 &
                -     gg(3,3) * gg(1,2) ** 2

           ggu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / detg
           ggu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / detg
           ggu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / detg
           ggu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / detg
           ggu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / detg
           ggu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / detg
           ggu(2,1) = ggu(1,2)
           ggu(3,1) = ggu(1,3)
           ggu(3,2) = ggu(2,3)

           gUP(1,1) = -one_over_alph_squared
           do a = 1, 3
              gUP(1,a + 1) = one_over_alph_squared * beta(a)
              gUP(a + 1,1) = one_over_alph_squared * beta(a)
           end do
           do a = 1, 3
              do b = 1, 3
                 gUP(a + 1, b + 1) = ggu(a,b) - one_over_alph_squared * beta(a) * beta(b)
              end do
           end do

           ! if (abs(x(i,j,k) - 0.8) <= 0.1 .and. abs(y(i,j,k) - 2.0) <= 0.1 .and. abs(z(i,j,k) - 3.0) <= 0.1) then
           !    write (*,*) "HERE detg", detg
           !    write (*,*) "HERE gg1", gg(1,:)
           !    write (*,*) "HERE gg2", gg(2,:)
           !    write (*,*) "HERE gg3", gg(3,:)
           !    write (*,*) "HERE ggu1", ggu(1,:)
           !    write (*,*) "HERE ggu2", ggu(2,:)
           !    write (*,*) "HERE ggu3", ggu(3,:)
           !    write (*,*) "HERE gab1", gab(1,:)
           !    write (*,*) "HERE gab2", gab(2,:)
           !    write (*,*) "HERE gab3", gab(3,:)
           !    write (*,*) "HERE gab4", gab(4,:)
           !    write (*,*) "HERE gUP1", gUP(1,:)
           !    write (*,*) "HERE gUP2", gUP(2,:)
           !    write (*,*) "HERE gUP3", gUP(3,:)
           !    write (*,*) "HERE gUP4", gUP(4,:)
           !    write (*,*) "HERE alph", alph
           !    write (*,*) "HERE dt_alph", dt_alph
           !    write (*,*) "HERE dalph", dalph(:)
           !    write (*,*) "HERE beta", beta(:)
           !    write (*,*) "HERE dbeta1", dbeta(1,:)
           !    write (*,*) "HERE dbeta2", dbeta(2,:)
           !    write (*,*) "HERE dbeta3", dbeta(3,:)
           !    write (*,*) "HERE dtbeta", dt_beta(:)
           !    write (*,*) "HERE dgg11", dgg(1,1,:)
           !    write (*,*) "HERE dgg12", dgg(1,2,:)
           !    write (*,*) "HERE dgg13", dgg(1,3,:)
           !    write (*,*) "HERE dgg21", dgg(2,1,:)
           !    write (*,*) "HERE dgg22", dgg(2,2,:)
           !    write (*,*) "HERE dgg23", dgg(2,3,:)
           !    write (*,*) "HERE dgg31", dgg(3,1,:)
           !    write (*,*) "HERE dgg32", dgg(3,2,:)
           !    write (*,*) "HERE dgg33", dgg(3,3,:)
           !    write (*,*) "HERE dt_gg1", dt_gg(1,:)
           !    write (*,*) "HERE dt_gg2", dt_gg(2,:)
           !    write (*,*) "HERE dt_gg3", dt_gg(3,:)
           !    write (*,*) "HERE dgab11", dgab(1,1,:)
           !    write (*,*) "HERE dgab12", dgab(1,2,:)
           !    write (*,*) "HERE dgab13", dgab(1,3,:)
           !    write (*,*) "HERE dgab14", dgab(1,4,:)
           !    write (*,*) "HERE dgab21", dgab(2,1,:)
           !    write (*,*) "HERE dgab22", dgab(2,2,:)
           !    write (*,*) "HERE dgab23", dgab(2,3,:)
           !    write (*,*) "HERE dgab24", dgab(2,4,:)
           !    write (*,*) "HERE dgab31", dgab(3,1,:)
           !    write (*,*) "HERE dgab32", dgab(3,2,:)
           !    write (*,*) "HERE dgab33", dgab(3,3,:)
           !    write (*,*) "HERE dgab34", dgab(3,4,:)
           !    write (*,*) "HERE dgab41", dgab(4,1,:)
           !    write (*,*) "HERE dgab42", dgab(4,2,:)
           !    write (*,*) "HERE dgab43", dgab(4,3,:)
           !    write (*,*) "HERE dgab44", dgab(4,4,:)

           !    stop
           ! end if

           ! Christoffel symbols with all indices downstairs
           cf1 = 0
           do a = 1, 4
             do b = 1, 4
               do c = b, 4
                 cf1(a,b,c) = 0.5d0 * (dgab(a,b,c) + dgab(a,c,b) - dgab(b,c,a))
               end do
             end do
           end do
           cf1(:,2,1) = cf1(:,1,2)
           cf1(:,3,1) = cf1(:,1,3)
           cf1(:,4,1) = cf1(:,1,4)
           cf1(:,3,2) = cf1(:,2,3)
           cf1(:,4,2) = cf1(:,2,4)
           cf1(:,4,3) = cf1(:,3,4)

           cf2 = 0
           do a = 1, 4
             do b = 1, 4
               do c = b, 4
                 do m = 1, 4
                   cf2(a,b,c) = cf2(a,b,c) + gUP(a,m) * cf1(m,b,c)
                 end do
               end do
             end do
           end do
           cf2(:,2,1) = cf2(:,1,2)
           cf2(:,3,1) = cf2(:,1,3)
           cf2(:,4,1) = cf2(:,1,4)
           cf2(:,3,2) = cf2(:,2,3)
           cf2(:,4,2) = cf2(:,2,4)
           cf2(:,4,3) = cf2(:,3,4)

           Gamma_ttt(i,j,k) = cf2(1,1,1)
           Gamma_ttx(i,j,k) = cf2(1,1,2)
           Gamma_tty(i,j,k) = cf2(1,1,3)
           Gamma_ttz(i,j,k) = cf2(1,1,4)
           Gamma_txx(i,j,k) = cf2(1,2,2)
           Gamma_txy(i,j,k) = cf2(1,2,3)
           Gamma_txz(i,j,k) = cf2(1,2,4)
           Gamma_tyy(i,j,k) = cf2(1,3,3)
           Gamma_tyz(i,j,k) = cf2(1,3,4)
           Gamma_tzz(i,j,k) = cf2(1,4,4)
           Gamma_xtt(i,j,k) = cf2(2,1,1)
           Gamma_xtx(i,j,k) = cf2(2,1,2)
           Gamma_xty(i,j,k) = cf2(2,1,3)
           Gamma_xtz(i,j,k) = cf2(2,1,4)
           Gamma_xxx(i,j,k) = cf2(2,2,2)
           Gamma_xxy(i,j,k) = cf2(2,2,3)
           Gamma_xxz(i,j,k) = cf2(2,2,4)
           Gamma_xyy(i,j,k) = cf2(2,3,3)
           Gamma_xyz(i,j,k) = cf2(2,3,4)
           Gamma_xzz(i,j,k) = cf2(2,4,4)
           Gamma_ytt(i,j,k) = cf2(3,1,1)
           Gamma_ytx(i,j,k) = cf2(3,1,2)
           Gamma_yty(i,j,k) = cf2(3,1,3)
           Gamma_ytz(i,j,k) = cf2(3,1,4)
           Gamma_yxx(i,j,k) = cf2(3,2,2)
           Gamma_yxy(i,j,k) = cf2(3,2,3)
           Gamma_yxz(i,j,k) = cf2(3,2,4)
           Gamma_yyy(i,j,k) = cf2(3,3,3)
           Gamma_yyz(i,j,k) = cf2(3,3,4)
           Gamma_yzz(i,j,k) = cf2(3,4,4)
           Gamma_ztt(i,j,k) = cf2(4,1,1)
           Gamma_ztx(i,j,k) = cf2(4,1,2)
           Gamma_zty(i,j,k) = cf2(4,1,3)
           Gamma_ztz(i,j,k) = cf2(4,1,4)
           Gamma_zxx(i,j,k) = cf2(4,2,2)
           Gamma_zxy(i,j,k) = cf2(4,2,3)
           Gamma_zxz(i,j,k) = cf2(4,2,4)
           Gamma_zyy(i,j,k) = cf2(4,3,3)
           Gamma_zyz(i,j,k) = cf2(4,3,4)
           Gamma_zzz(i,j,k) = cf2(4,4,4)

           if (save_dgab/= 0) then
              ! Here the derivative index is the first, but in dgab is the last
              ! (We follow the same names as Gammas)
              g_ttt(i,j,k) = dgab(1,1,1)
              g_ttx(i,j,k) = dgab(1,2,1)
              g_tty(i,j,k) = dgab(1,3,1)
              g_ttz(i,j,k) = dgab(1,4,1)
              g_txx(i,j,k) = dgab(2,2,1)
              g_txy(i,j,k) = dgab(2,3,1)
              g_txz(i,j,k) = dgab(2,4,1)
              g_tyy(i,j,k) = dgab(3,3,1)
              g_tyz(i,j,k) = dgab(3,4,1)
              g_tzz(i,j,k) = dgab(4,4,1)
              g_xtt(i,j,k) = dgab(1,1,2)
              g_xtx(i,j,k) = dgab(1,2,2)
              g_xty(i,j,k) = dgab(1,3,2)
              g_xtz(i,j,k) = dgab(1,4,2)
              g_xxx(i,j,k) = dgab(2,2,2)
              g_xxy(i,j,k) = dgab(2,3,2)
              g_xxz(i,j,k) = dgab(2,4,2)
              g_xyy(i,j,k) = dgab(3,3,2)
              g_xyz(i,j,k) = dgab(3,4,2)
              g_xzz(i,j,k) = dgab(4,4,2)
              g_ytt(i,j,k) = dgab(1,1,3)
              g_ytx(i,j,k) = dgab(1,2,3)
              g_yty(i,j,k) = dgab(1,3,3)
              g_ytz(i,j,k) = dgab(1,4,3)
              g_yxx(i,j,k) = dgab(2,2,3)
              g_yxy(i,j,k) = dgab(2,3,3)
              g_yxz(i,j,k) = dgab(2,4,3)
              g_yyy(i,j,k) = dgab(3,3,3)
              g_yyz(i,j,k) = dgab(3,4,3)
              g_yzz(i,j,k) = dgab(4,4,3)
              g_ztt(i,j,k) = dgab(1,1,4)
              g_ztx(i,j,k) = dgab(1,2,4)
              g_zty(i,j,k) = dgab(1,3,4)
              g_ztz(i,j,k) = dgab(1,4,4)
              g_zxx(i,j,k) = dgab(2,2,4)
              g_zxy(i,j,k) = dgab(2,3,4)
              g_zxz(i,j,k) = dgab(2,4,4)
              g_zyy(i,j,k) = dgab(3,3,4)
              g_zyz(i,j,k) = dgab(3,4,4)
              g_zzz(i,j,k) = dgab(4,4,4)
           end if
        end do
     end do
  end do

9999 continue

end subroutine Christoffel_compute
