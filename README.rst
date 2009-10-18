===============================================================================
cysolar 
===============================================================================

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
port of pysolar_ to cython_
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

all the brain work was done by Brandon of pysolar_ this is a naive port
for speed. Thanks to Brandon for allowing this to be released under A 
BSD license.

example usage ::

    >>> from cysolar import get_radiation
    >>> import datetime
    >>> d = datetime.datetime(2007, 2, 18, 20, 13, 1, 130320)
    >>> lon, lat = -71.382, 42.206
    >>> get_radiation(d, lon, lat)
    803.57188602589667



.. _cython: http://cython.org/
.. _pysolar: http://github.com/pingswept/pysolar
