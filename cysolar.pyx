"""
cysolar
=======

port of pysolar to cython. see README.rst

(License is BSD, thanks to Brandon)

"""
DEF DELTA_SECONDS = 65.0
DEF PI = 3.1415926535897931
DEF EARTH_RADIUS = 6378140.0 # meters

cdef extern from "solar_constants.h":
    double L0[64][3]
    double L1[34][3]
    double L2[20][3]
    double L3[7][3]
    double L4[3][3]
    double L5[1][3]

    double B0[5][3]
    double B1[2][3]

    double R0[40][3]
    double R1[10][3]
    double R2[6][3]
    double R3[2][3]
    double R4[1][3]

    double nutation_coefficients[63][4]
    int aberration_sin_terms[63][5]

    double aberration_poly(double jce, double a, double b, double c, double d)

    int get_day_of_year(int y, int m, int d)
    int months[12]

cdef extern from "math.h":
    double sin(double x)
    double exp(double x)
    double atan(double x)
    double tan(double x)
    double cos(double x)
    double fmod(double x, double y)
    double atan2(double x, double y)
    double asin(double x)
    double floor(double x)
    double abs(double x)

cdef inline double get_julian_century(double julian_day):
    return (julian_day - 2451545.0) / 36525.0

cdef inline double get_julian_day(int year, unsigned int month, 
                                  unsigned int day, unsigned int hour, 
                                  unsigned int minute, unsigned int second, 
                                  unsigned int msecond):   
    # based on NREL/TP-560-34302 by Andreas and Reda. does not
    # accept years before 0 because of bounds check on Python's datetime.year field
    cdef double julian_day, dday, gregorian_offset
    if month <= 2.0 :       # shift to accomodate leap years?
        year -=  1
        month += 12

    dday = day + (((hour * 3600.0) + (minute * 60.0) +\
                   second + (msecond / 1000000.0)) / 86400.0)

    julian_day = floor(365.25 * (year + 4716.0)) + \
            floor(30.6001 * (month + 1.0)) + dday - 1524.5
    if julian_day <= 2299160.0:
        return julian_day # before October 5, 1852

    gregorian_offset = 2.0 - (year // 100.0) + ((year // 100.0) / 4.0)
    return julian_day + gregorian_offset # after October 5, 1852


cdef inline double get_julian_ephemeris_century(double julian_ephemeris_day):
    return (julian_ephemeris_day - 2451545.0) / 36525.0

cdef inline double get_julian_ephemeris_day(double julian_day):
    """delta_seconds is value referred to by astronomers as Delta-T, defined as the difference between
    Dynamical Time (TD) and Universal Time (UT). In 2007, it's around 65 seconds.
    A list of values for Delta-T can be found here: ftp://maia.usno.navy.mil/ser7/deltat.data"""
    return julian_day + (DELTA_SECONDS / 86400.0)

cdef inline double get_julian_ephemeris_millenium(double julian_ephemeris_century):
    return (julian_ephemeris_century / 10.0)


cdef inline double radians(double d):
    return PI * d / 180.0

cdef inline double degrees(double r):
    return 180.0 / PI * r

cdef inline double air_mass_ratio(double altitude_deg):
    # from Masters, p. 412
    # warning: pukes on input of zero
    return 1.0 / sin(radians(altitude_deg))

cdef inline double apparent_extraterrestrial_flux(unsigned int day):
    # from Masters, p. 412 W / m*m
    return 1160.0 + (75.0 * sin((360.0/365.0) * (day - 275.0)))

# Renewable and efficient electric power systems By Gilbert M. Masters
cdef inline double optical_depth(unsigned int day):
    # from Masters, p. 412
    return 0.174 + (0.035 * sin((360.0 / 365) * (day - 100.0)))

cdef inline double get_flattened_latitude(double latitude):
    cdef double latitude_rad = radians(latitude)
    return degrees(atan(0.99664719 * tan(latitude_rad)))

cdef double get_projected_radial_distance(double elevation, double latitude):
    cdef double flattened_latitude_rad = radians(get_flattened_latitude(latitude))
    cdef double latitude_rad = radians(latitude)
    return cos(flattened_latitude_rad) + (elevation * cos(latitude_rad) / EARTH_RADIUS)

cdef double get_projected_axial_distance(double elevation, double latitude):
    cdef double flattened_latitude_rad = radians(get_flattened_latitude(latitude))
    cdef double latitude_rad = radians(latitude)
    return 0.99664719 * sin(flattened_latitude_rad) + (elevation * sin(latitude_rad) / EARTH_RADIUS)

cdef inline double get_aberration_correction(double radius_vector):    
    # r is earth radius vector [astronomical units]
   return -20.4898/(3600.0 * radius_vector)

cdef inline double get_equatorial_horizontal_parallax(double radius_vector):
    return 8.794 / (3600.0 / radius_vector)

cdef double get_nutation_aberration_XY(double jce, unsigned int i, double *x):
    cdef int nx = 5 # hard code, but this is from the arr in aberration_sin_terms
    cdef double sigmaxy = 0.0
    cdef unsigned int j
    cdef int *row = <int * >aberration_sin_terms
    for j in range(nx):
        sigmaxy += x[j] * row[nx * i + j]
    return sigmaxy

cdef double* precalculate_aberrations(double jce):
    cdef double x[5]
    # test against: http://forums.swishzone.com/lofiversion/index.php?t7814-50.html
    x[0] = aberration_poly(jce, 297.85036, 445267.111480,  -0.0019142, 189474.0)
    x[1] = aberration_poly(jce, 357.52772, 35999.050340, -0.0001603, -300000.0)
    x[2] = aberration_poly(jce, 134.96298, 477198.867398, 0.0086972, 56250.0)
    x[3] = aberration_poly(jce, 93.27191, 483202.017538, -0.0036825, 327270.0)
    x[4] = aberration_poly(jce, 125.04452, -1934.136261, 0.0020708, 450000.0)
    return x

cdef double* get_nutations(double jde):
    cdef unsigned int i, nc = 63
    cdef double jce = get_julian_ephemeris_century(jde)
    cdef double nutation_long = 0, nutation_oblique = 0, sigmaxy

    cdef double *x = precalculate_aberrations(jce)
    cdef double nutations[2]

    cdef double *nc1d = <double *>nutation_coefficients
    for i in range(nc):

        sigmaxy = radians(get_nutation_aberration_XY(jce, i, x))
        nutation_long += (nc1d[i* 4 + 0] \
                       + (nc1d[i * 4 + 1] * jce)) * sin(sigmaxy)

        nutation_oblique += (nc1d[i * 4 + 2] \
                          + (nc1d[i * 4 + 3] * jce)) * cos(sigmaxy)
    # 36000000 scales from 0.0001 arcseconds to degrees
    nutations[0] = nutation_long / 36000000.0
    nutations[1] = nutation_oblique / 36000000.0

    return nutations

cdef double get_coefficient(double jme, double **constant_array, unsigned int alen):
    cdef double c = 0.0
    cdef unsigned int i
    # hopefully i'm doing this pointer arith correctly...
    cdef double *a1d = <double *>constant_array
    
    for i in range(alen):
        c += a1d[3 * i + 0] * cos(a1d[3 * i + 1]  + a1d[3 * i + 2] * jme) 
    return c

cdef inline double get_heliocentric_latitude(double jme):
    cdef double b0 = get_coefficient(jme, <double **>B0, 5)
    cdef double b1 = get_coefficient(jme, <double **>B1, 2)
    return degrees((b0 + (b1 * jme)) / 10 ** 8)

cdef double get_heliocentric_longitude(double jme):
    cdef double l0 = get_coefficient(jme, <double **>L0, 64)
    cdef double l1 = get_coefficient(jme, <double **>L1, 34)
    cdef double l2 = get_coefficient(jme, <double **>L2, 20)
    cdef double l3 = get_coefficient(jme, <double **>L3, 7)
    cdef double l4 = get_coefficient(jme, <double **>L4, 3)
    cdef double l5 = get_coefficient(jme, <double **>L5, 1)

    cdef double l = (l0 + l1 * jme + l2 * jme ** 2 + l3 * jme ** 3 + l4 * jme ** 4 + l5 * jme ** 5) / 10 ** 8
    return degrees(l) % 360

cdef double get_radius_vector(double jme):
   cdef double r0 = get_coefficient(jme, <double **>R0, 40)
   cdef double r1 = get_coefficient(jme, <double **>R1, 10)
   cdef double r2 = get_coefficient(jme, <double **>R2, 6)
   cdef double r3 = get_coefficient(jme, <double **>R3, 2)
   cdef double r4 = get_coefficient(jme, <double **>R4, 1)

   return (r0 + r1 * jme + r2 * jme ** 2 + r3 * jme ** 3 + r4 * jme ** 4) / 10 ** 8

cdef inline double get_mean_sidereal_time(double julian_day):
   # This function doesn't agree with Andreas and Reda as well as it should. Works to ~5 sig figs in current unit test
   cdef double jc = get_julian_century(julian_day)
   cdef double sidereal_time =  280.46061837 + (360.98564736629 * (julian_day - 2451545.0)) + (0.000387933 * jc ** 2) - (jc ** 3 / 38710000)
   return fmod(sidereal_time, 360.0)

cdef double get_true_ecliptic_obliquity(double jme, double *nutation):
   cdef double u = jme/10.0
   cdef double mean_obliquity = 84381.448 - (4680.93 * u) - (1.55 * u ** 2) +\
           (1999.25 * u ** 3) - (51.38 * u ** 4) -(249.67 * u ** 5) -\
           (39.05 * u ** 6) + (7.12 * u ** 7) + (27.87 * u ** 8) +\
           (5.79 * u ** 9) + (2.45 * u ** 10)
   return (mean_obliquity / 3600.0) + nutation[1]

cdef double get_apparent_sidereal_time(double julian_day, double jme, double *nutation, double true_eliptic_obliquity):
    # nutation is 2array of long, oblique
    return get_mean_sidereal_time(julian_day) + nutation[0] * \
            cos(true_eliptic_obliquity)

cdef inline double get_geocentric_latitude(double jme):
   return -1 * get_heliocentric_latitude(jme)

cdef inline double get_geocentric_longitude(double jme):
   return fmod((get_heliocentric_longitude(jme) + 180), 360.0)

cdef inline double get_apparent_sun_longitude(double geocentric_longitude, double *nutations, double ab_correction):
    return geocentric_longitude + nutations[0] + ab_correction

cdef double get_geocentric_sun_right_ascension(double apparent_sun_longitude, 
                                               double true_ecliptic_obliquity, 
                                               double geocentric_latitude):
    cdef double apparent_sun_longitude_rad = radians(apparent_sun_longitude)
    cdef double true_ecliptic_obliquity_rad = radians(true_ecliptic_obliquity)
    cdef double geocentric_latitude_rad = radians(geocentric_latitude)

    cdef double a = sin(apparent_sun_longitude_rad) * cos(true_ecliptic_obliquity_rad)
    cdef double b = tan(geocentric_latitude_rad) * sin(true_ecliptic_obliquity_rad)
    cdef double c = cos(apparent_sun_longitude_rad)
    cdef double alpha = atan2((a - b),  c)


    cdef double rem = fmod(degrees(alpha), 360.0)
    # the python version doesnt have to do this, CHECK!!!
    return rem if rem > 0 else rem + 360


cdef inline double get_local_hour_angle(double apparent_sidereal_time, double longitude, double geocentric_sun_right_ascension):
    cdef double rem = fmod((apparent_sidereal_time + longitude -\
                 geocentric_sun_right_ascension), 360.0)
    # the python version doesnt have to do this, CHECK!!!
    return rem if rem > 0 else rem + 360

cdef double get_parallax_sun_right_ascension(double prd, 
                                             double equatorial_horizontal_parallax, 
                                             double local_hour_angle, 
                                             double geocentric_sun_declination):
    cdef double ehp_rad = radians(equatorial_horizontal_parallax)
    cdef double lha_rad = radians(local_hour_angle)
    cdef double gsd_rad = radians(geocentric_sun_declination)
    cdef double a = -1 * prd * sin(ehp_rad) * sin(lha_rad)
    cdef double b =  cos(gsd_rad) - prd * sin(ehp_rad) * cos(lha_rad)
    cdef double parallax = atan2(a, b)
    return degrees(parallax)

cdef inline double get_topocentric_local_hour_angle(double local_hour_angle, 
                                                    double parallax_sun_right_ascension):
    return local_hour_angle - parallax_sun_right_ascension

cdef double get_topocentric_sun_declination(double geocentric_sun_declination, 
                                            double projected_axial_distance, 
                                            double equatorial_horizontal_parallax, 
                                            double parallax_sun_right_ascension, 
                                            double local_hour_angle):
    cdef double gsd_rad = radians(geocentric_sun_declination)
    cdef double pad = projected_axial_distance
    cdef double ehp_rad = radians(equatorial_horizontal_parallax)
    cdef double psra_rad = radians(parallax_sun_right_ascension)
    cdef double lha_rad = radians(local_hour_angle)
    cdef double a = (sin(gsd_rad) - pad * sin(ehp_rad)) * cos(psra_rad)
    cdef double b = cos(gsd_rad) - (pad * sin(ehp_rad) * cos(lha_rad))
    return degrees(atan2(a, b))

cdef double get_geocentric_sun_declination(double apparent_sun_longitude, double true_ecliptic_obliquity, double geocentric_latitude):
    cdef double apparent_sun_longitude_rad = radians(apparent_sun_longitude)
    cdef double true_ecliptic_obliquity_rad = radians(true_ecliptic_obliquity)
    cdef double geocentric_latitude_rad = radians(geocentric_latitude)

    cdef double a = sin(geocentric_latitude_rad) * cos(true_ecliptic_obliquity_rad)
    cdef double b = cos(geocentric_latitude_rad) * sin(true_ecliptic_obliquity_rad) * sin(apparent_sun_longitude_rad)
    cdef double delta = asin(a + b)
    return degrees(delta)

cdef double get_topocentric_elevation_angle(double latitude, 
                                            double topocentric_sun_declination, 
                                            double topocentric_local_hour_angle):
    cdef double latitude_rad = radians(latitude)
    cdef double tsd_rad = radians(topocentric_sun_declination)
    cdef double tlha_rad = radians(topocentric_local_hour_angle)
    return degrees(asin((sin(latitude_rad) * sin(tsd_rad)) + cos(latitude_rad) * cos(tsd_rad) * cos(tlha_rad)))

cdef inline double get_refraction_correction(double pressure_millibars, 
                                             double temperature_celsius,
                                             double topocentric_elevation_angle):
    cdef double tea = topocentric_elevation_angle
    cdef double temperature_kelvin = temperature_celsius + 273.15
    cdef double a = pressure_millibars * 283.0 * 1.02
    cdef double b = 1010.0 * temperature_kelvin * 60.0 *\
                        tan(radians(tea + (10.3/(tea + 5.11))))
    return a / b



# NOTE: this is lon, lat order where as the python is lat, lon
# to maintain compatibility with pysolar.
cdef double cget_altitude(double lon, double lat, unsigned int year,
                           unsigned int month,
                           unsigned int day,
                           unsigned int hour,
                           unsigned int minute,
                           unsigned int second,
                           unsigned int msecond,
                           double elevation,
                           double temperature_celsius,
                         double pressure_millibars):
    cdef double jd, jde, jce, jme
    cdef double geocentric_latitude, geocentric_longitude, radius_vector, aberration_correction
    cdef double projected_radial_distance, projected_axial_distance, equatorial_horizontal_parallax
    cdef double *nutations
    cdef double flattened_latitude_rad = radians(get_flattened_latitude(lat))
    cdef double latitude_rad = radians(lat)
    cdef double apparent_sidereal_time, true_eliptic_obliquity
    cdef double apparent_sun_longitude, geocentric_sun_right_ascension
    cdef double geocentric_sun_declination, local_hour_angle, parallax_sun_right_ascension
    cdef double topocentric_local_hour_angle, topocentric_sun_declination
    cdef double topocentric_elevation_angle, refraction_correction, true_ecliptic_obliquity
    projected_radial_distance = cos(flattened_latitude_rad) + (elevation * cos(latitude_rad) / EARTH_RADIUS)
    projected_axial_distance = 0.99664719 * sin(flattened_latitude_rad) + (elevation * sin(latitude_rad) / EARTH_RADIUS)

    jd = get_julian_day(year, month, day, hour, minute, second, msecond)
    jde = get_julian_ephemeris_day(jd) #, 65)
    jce = get_julian_ephemeris_century(jde)
    jme = get_julian_ephemeris_millenium(jce)

    geocentric_longitude = get_geocentric_longitude(jme)
    geocentric_latitude = get_geocentric_latitude(jme)

    radius_vector = get_radius_vector(jme)
    aberration_correction = get_aberration_correction(radius_vector)
    equatorial_horizontal_parallax = get_equatorial_horizontal_parallax(radius_vector)

    # long, oblique nutations
    nutations = get_nutations(jde)
    true_ecliptic_obliquity = get_true_ecliptic_obliquity(jme, nutations)
    apparent_sidereal_time = get_apparent_sidereal_time(jd, jme, nutations, true_ecliptic_obliquity)

    # calculations dependent on location and time
    apparent_sun_longitude = get_apparent_sun_longitude(geocentric_longitude,
                                                        nutations, 
                                                        aberration_correction)

    geocentric_sun_right_ascension = get_geocentric_sun_right_ascension(
                                                     apparent_sun_longitude,
                                                     true_ecliptic_obliquity,
                                                          geocentric_latitude)

    geocentric_sun_declination = get_geocentric_sun_declination(apparent_sun_longitude, true_ecliptic_obliquity, geocentric_latitude) 
    local_hour_angle = get_local_hour_angle(apparent_sidereal_time, lon, geocentric_sun_right_ascension)


    parallax_sun_right_ascension = get_parallax_sun_right_ascension(projected_radial_distance, equatorial_horizontal_parallax, local_hour_angle, geocentric_sun_declination)
    topocentric_local_hour_angle = get_topocentric_local_hour_angle(local_hour_angle, parallax_sun_right_ascension)
    topocentric_sun_declination = get_topocentric_sun_declination(geocentric_sun_declination, projected_axial_distance, equatorial_horizontal_parallax, parallax_sun_right_ascension, local_hour_angle)

    topocentric_elevation_angle = get_topocentric_elevation_angle(lat, topocentric_sun_declination, topocentric_local_hour_angle)

    refraction_correction = get_refraction_correction(pressure_millibars, temperature_celsius, topocentric_elevation_angle)

    return topocentric_elevation_angle + refraction_correction

cdef inline double get_air_mass_ratio(double altitude_deg):
    # from Masters, p. 412
    # warning: pukes on input of zero
    return 1.0 / sin(radians(altitude_deg))

cdef inline double get_apparent_extraterrestrial_flux(double day):
    # from Masters, p. 412
    return 1160 + (75 * sin((360.0/365.0) * (day - 275)))

cdef inline double get_optical_depth(double day):
    # from Masters, p. 412
    return 0.174 + (0.035 * sin((360.0/365.0) * (day - 100)))

import datetime
def py_get_day_of_year(utc_datetime):
    year_start = datetime.datetime(utc_datetime.year, 1, 1)
    delta = (utc_datetime - year_start)
    return delta.days


def get_altitude(double lat, double lon, 
                 object adatetime, 
                 double elevation=0.0,
                 double temperature_celsius=25, 
                 double pressure_millibars=1013.25):
    d = adatetime
    return cget_altitude(lon, lat, d.year, d.month, d.day, d.hour,
                         d.minute, d.second, d.microsecond, elevation,
                         temperature_celsius, pressure_millibars)

def get_radiation(adatetime, double lon, double lat, 
                  double elevation=0.0,
                  double temperature_celsius=25,
                  double pressure_millibars=1013.25):
    """
    given a date, lon and lat, return the radiation in watts/m**2
    :param adatetime: a python datetime
    :param lon: longitude
    :param lat: latitude
    :param elevation: the elevation in meters.
    :param temperature_celsius: ...
    :param pressure_millibars: ...
    """
    d = adatetime
    cdef double alt_deg = cget_altitude(lon, lat, d.year, d.month,
                                        d.day, d.hour, d.minute, d.second,
                                        d.microsecond, elevation,
                                        temperature_celsius,
                                        pressure_millibars)
    cdef unsigned int day = py_get_day_of_year(adatetime)
    return _radiation(day, alt_deg)


def get_radiation_for_year(year, lon, lat, elevation):
    cdef double temp_c = 25.0
    cdef double pressure = 1013.0
    cdef unsigned month, day, day_of_y = 0, hour, minute, i

    cdef double rads[366]
    cdef double daily_rad, alt_deg, flux, optical_depth, air_mass_ratio

    ops = []
    fluxs = []
    for month in range(12):
        for day in range(1, months[month] + 1):
            day_of_y += 1
            daily_rad = 0.0;
            flux = get_apparent_extraterrestrial_flux(day_of_y)
            optical_depth = get_optical_depth(day_of_y)
            ops.append(optical_depth)
            fluxs.append(flux)
            print "%.3f, %.3f" % (flux, optical_depth)
            for hour in range(0, 24):
                for minute in range(0, 60):
                    alt_deg = cget_altitude(lon, lat, year, month,
                                        day, hour, minute, 0,
                                        0, elevation,
                                        temp_c,
                                        pressure)
                    if alt_deg > 0: 
                        air_mass_ratio = get_air_mass_ratio(alt_deg)
                        daily_rad += (flux * exp(-1.0 * optical_depth * \
                                             air_mass_ratio))
            rads[day_of_y - 1] = daily_rad / (24 * 60)
    return [rads[i] for i in range(365)], ops, fluxs

cdef double _radiation(unsigned int day, double altitude_deg):
    if altitude_deg <= 0: return 0.0
    cdef double flux = get_apparent_extraterrestrial_flux(day)
    cdef double optical_depth = get_optical_depth(day)
    cdef double air_mass_ratio = get_air_mass_ratio(altitude_deg)
    return flux * exp(-1.0 * optical_depth * air_mass_ratio)

cpdef get_radiation_direct(utc_datetime, double altitude_deg):
    """
    see `get_radiation` for simpler use.
    takes a python datetime object and an altitude in degrees
    from `get_altitude` and returns the radiation in Watts / meter**2
    see: http://wiki.github.com/pingswept/pysolar/examples

    >>> import datetime
    >>> from cysolar import get_radiation_direct, get_altitude
    >>> d = datetime.datetime(2007, 2, 18, 20, 13, 1, 130320)
    >>> lon, lat = -71.382, 42.206
    >>> altitude_deg = get_altitude(lat, lon, d)
    >>> altitude_deg = get_altitude(42.206, -71.382, d)
    >>> altitude_deg
    20.373453131123568

    >>> get_radiation_direct(d, altitude_deg)
    803.57188602589667

    """

    # from Masters, p. 412
    cdef unsigned int day
    if altitude_deg <= 0: return 0.0
    day = py_get_day_of_year(utc_datetime)
    return _radiation(day, altitude_deg)
