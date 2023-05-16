import arcpy
from arcpy import sa
from arcpy.sa import *
from arcpy import env

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True


# GPP = GPP_shade + GPP_sun

for yr in range (1981, 2019):
    forest_map = arcpy.Raster("")
    forest_tmp = Con(IsNull(forest_map), 0, forest_map)
    # eps shade
    forest_eps_shade = Con(forest_tmp == 1, 3.75,
                           Con(forest_tmp == 2, 3.26, Con(forest_tmp == 3, 3.40, Con(forest_tmp == 4, 3.00, 0))))
    forest_eps_shade = Con(forest_tmp == 1, 0.92,
                           Con(forest_tmp == 2, 1.44, Con(forest_tmp == 3, 0.89, Con(forest_tmp == 4, 0.80, 0))))
    forest_Topt = Con(forest_tmp == 1, 23.1,
                           Con(forest_tmp == 2, 25.8, Con(forest_tmp == 3, 19.7, Con(forest_tmp == 4, 24.5, 0))))
    forest_albedo = Con(forest_tmp == 1, 0.18,
                      Con(forest_tmp == 2, 0.18, Con(forest_tmp == 3, 0.15, Con(forest_tmp == 4, 0.17, 0))))
    forest_omega = Con(forest_tmp == 1, 0.8,
                        Con(forest_tmp == 2, 0.8, Con(forest_tmp == 3, 0.6, Con(forest_tmp == 4, 0.7, 0))))

    for d in range(1, 47):
        # Calculate LAIshape and LAIsun
        # input solar zenith angle; clumping index; LAI
        LAI = arcpy.Raster("")
        omega = 0.0
        costha_m = 0.0

        LAIsun = 2 * costha_m * (1 - Exp(-0.5 * omega * LAI / costha_m))
        LAIshade = LAI - LAIsun

        # Calculate PARdir, Pdif
        PAR = arcpy.Raster("")
        R = PAR / (1367 * 0.43 * costha_m)
        PARdif = PAR * (0.7527 + 3.8453 * R - 16.316 * R * R + 18.962 * R * R * R - 7.0802 * R * R * R * R)
        PARdir = PAR - PARdif

        cos_szaa = 0.537 + 0.025 * LAI
        PARdif_under = PARdif * Exp(-0.5 * omega * LAI / cos_szaa)
        C = 0.07 * omega * PARdir * (1.1 - 0.1 * LAI) * Exp(-costha_m)

        # Calculate APAR
        alpha = 0.15

        PARshade = (1 - alpha) * ((PARdif - PARdif_under) / LAI + C)
        PARsun = (1 - alpha) * PARdir * 0.5 / costha_m + PARshade
        APARsun = LAIsun * PARsun
        APARshade = LAIshade * PARshade

        #  Temperature scalar
        tmp_avg = arcpy.Raster("")

        Tmin = 0
        Tmax = 40
        Topt = 0.0

        Ts_tmp = (tmp_avg - Tmin) * (tmp_avg - Tmax) / (
                    (tmp_avg - Tmin) * (tmp_avg - Tmax) - (tmp_avg - Topt) * (tmp_avg - Topt))
        Ts = Con(tmp_avg < 0, 0, Con(tmp_avg > 40, 0, Ts_tmp))

        # VPD scalar
        vpd = arcpy.Raster("")
        vpd_min = 930
        vpd_max = 4100

        vpds_tmp = (vpd_max - vpd) / (vpd_max - vpd_min)
        vpds = Con(vpd < vpd_min, 1, Con(vpd > vpd_max, 0, vpds_tmp))

        # CO2 effects
        # input temperature vpd Ca
        Rgas = 8.314
        ita = 0.8903
        P = 21
        Ca = 0
        fi = 4.22 * Exp((37830 * (tmp_avg - 298.15)) / (298.15 * tmp_avg * 8.314))
        Kc = 39.97 * Exp(79.43 * (tmp_avg - 298.15) / (298.15 * Rgas * tmp_avg))
        Ko = 27480 * Exp(36.38 * (tmp_avg - 298.15) / (298.15 * Rgas * tmp_avg))
        K = Kc * (1 + P / Ko)
        kesi = SquareRoot(356.51 * K / (1.6 * ita))
        gamma = kesi / (kesi + SquareRoot(vpd))
        Ci = Ca * gamma
        CO2s = (Ci - fi) / (Ci + 2 * fi)

        #
        eps_sun = 0.0
        eps_shade = 0.0

        GPP_sun = eps_sun * APARsun * Ts * vpds * CO2s
        GPP_shade = eps_shade * APARshade * Ts * vpds * CO2s
        GPP = GPP_sun + GPP_shade
        GPP.save("")
