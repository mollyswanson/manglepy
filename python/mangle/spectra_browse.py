#!/usr/bin/env python
"""
Open a sequence of spectra in the new BOSS spectra browser.
"""
#Initial Version: jkp 2010.10.01

import webbrowser

url='https://spectra.sdss3.org:8000/spectrumDetail?plateid=PLATE&mjd=MJD&fiber=FIBER&reduction2d=v5_4_14'

def open_spectra(data):
    browser = webbrowser.get()
    print 'Plate Fiber MJD'
    for x in data:
        print x['PLATE'],x['MJD'],x['FIBERID']
        temp = url.replace('FIBER',str(x['FIBERID'])).replace('MJD',str(x['MJD'])).replace('PLATE',str(x['PLATE']))
        # if the browser is Safari, new=1 keeps things in the same window.
        browser.open(temp)
        raw_input('Press enter for next plot, or CTL-C to quit...')
#...
