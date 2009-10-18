import datetime
import sys
sys.path.insert(0, "../pysolar")
import solar
GetRadiationDirect = solar.radiation.GetRadiationDirect
GetAltitude = solar.GetAltitude

from cysolar import get_radiation_direct, get_altitude, get_radiation
d = datetime.datetime(2007, 2, 18, 20, 13, 1, 130320)
lon, lat = -71.382, 42.206
altitude_deg = get_altitude(lat, lon, d)
altitude_deg = get_altitude(lat, lon, d)
print altitude_deg
print GetAltitude(lat, lon, d)

print get_radiation_direct(d, altitude_deg)
print GetRadiationDirect(d, altitude_deg)

print get_radiation(d, lon, lat)
