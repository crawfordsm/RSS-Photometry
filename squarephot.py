
import numpy as np
from astropy.stats import sigma_clip

from pylab import *

def square_phot(data, xc, yc, xw, yw):
    """Return the number of counts from within an aperture.  The position of the object
       and the radius of interest should be provided
    
       Parameters
       ----------
       data: ~numpy.ndarray
           Array of counts

       xc:  float
           x-position for the source

       yc:  float
           x-position for the source

       xw:  float

           x-width for the source

       yw:  float
           y-width for the source


       Returns
       -------
       counts: float
           Number of counts within an aperture of size radius.  No
           corrections are applied for sub-pixel sampling.
    """
    y,x = np.indices(data.shape)
    mask = (x > xc-xw*0.5) * (x < xc+xw*0.5) * (y > yc-yw*0.5) * (y < yc+yw*0.5)
    return data[mask].sum()


if __name__=='__main__':
   import sys
   from astropy.io import fits

   if len(sys.argv) < 4:
      print """\nTo run program from the command line:
  
python aperturephot.py [img] [xc] [yc] [xw] [yw]

img -- name of image to run the program
xc  -- xcenter for source
yc  -- ycenter for source
radius -- radius to measure the total number of counts on
"""
      exit(1)

   hdu = fits.open(sys.argv[1])
   xc = float(sys.argv[2])    
   yc = float(sys.argv[3])    
   xw = float(sys.argv[4])    
   yw = float(sys.argv[5])    

   #remove the background

   #determine the photmetry
   counts = square_phot(hdu[0].data, xc, yc, xw, yw)

   #report the results
   print "No background Subtraction: ", counts
