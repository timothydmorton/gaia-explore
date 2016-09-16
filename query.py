#!/usr/bin/env python

from __future__ import print_function, division

import os, os.path

import pandas as pd
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u

catalogs = ['2mass', 'allwise', 'sdss9']

data = pd.read_hdf('data/merged_matches.h5','df')

i=1
for _, s in data.iterrows():
    name = s.epic_name
    print('{} ({} of {}).'.format(name, i, len(data)))
    folder = os.path.join('starmodels', name)
    if not os.path.exists(folder):
        os.makedirs(folder)
    ra,dec = s.epic_ra, s.epic_dec
    q = Vizier.query_region(SkyCoord(ra,dec, unit=(u.deg, u.deg)), 
                    radius=5*u.arcsec, catalog=catalogs)
    for t in q:
        t.write(os.path.join(folder,'{}.csv'.format(t._meta['ID'])), format='ascii.csv')
    i += 1
