
import glob
import numpy as np
from astropy.stats import sigma_clip

from aperturephot import aperture_phot
from pylab import *


if __name__=='__main__':
   import sys
   from astropy.io import fits


   outfile = sys.argv[1]
   radius = float(sys.argv[2])


   fout= open(outfile, 'w')

   #grab all the files
   infiles = glob.glob('*fit')

   #remove the background
   counts = []
   for img in infiles[0:]:
       hdu = fits.open(img)

       #background subtract the data
       bkgrd = np.median(sigma_clip(hdu[0].data, sig=3, iters=5, cenfunc=np.median, varfunc=np.var)[0])
       data = hdu[0].data - bkgrd

       #find the maximum pixel
       i = data.argmax() 
       ys, xs = data.shape
       yc = int(i/ys)
       xc = i%ys
      
       #perform photometry around that
       c = aperture_phot(hdu[0].data-bkgrd, xc, yc, radius)
       fout.write('%s %s %i %i %f\n' % (img, hdu[0].header["DATE-OBS"], xc, yc, c))
       counts.append(c)

   #plot the measurements
   plot(counts)
   show() 
      
