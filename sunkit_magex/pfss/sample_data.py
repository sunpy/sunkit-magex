"""
Functions to get sample data. The data is fetched and cached automatically
using `sunpy.data.manager`.
"""

from sunpy.data import manager


@manager.require('gong_map',
                 'https://gong2.nso.edu/oQR/zqs/202009/mrzqs200901/mrzqs200901t1304c2234_022.fits.gz',
                 'aad927d8f617f32b72255b862c4910f13640fc7ca13edf98288cd0735a2db6a0')
def get_gong_map():
    """
    Automatically download and unzip a sample GONG synoptic map.
    """
    return manager.get('gong_map')


@manager.require('adapt_map',
                 'https://gong.nso.edu/adapt/maps/gong/2020/adapt40311_044012_202010100000_i00005600n1.fts.gz',
                 'b8a04ec57eb3e797fe3de88f0d2b58b1a93607a914625341bb5a79268cc6314a')
def get_adapt_map():
    return manager.get('adapt_map')
