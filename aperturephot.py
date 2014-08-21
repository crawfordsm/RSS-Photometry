
import numpy as np
from astropy.stats import sigma_clip

from pylab import *

def aperture_phot(data, xc, yc, radius):
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

       radius:  float
           Within this radius, all sources will be summed

       Returns
       -------
       counts: float
           Number of counts within an aperture of size radius.  No
           corrections are applied for sub-pixel sampling.
    """
    y,x = np.indices(data.shape)
    r = ((x-xc)**2 + (y-yc)**2)**0.5
    return data[r <= radius].sum()


if __name__=='__main__':
   import sys
   from astropy.io import fits

   if len(sys.argv) < 4:
      print """\nTo run program from the command line:
  
python aperturephot.py [img] [xc] [yc] [radius]

img -- name of image to run the program
xc  -- xcenter for source
yc  -- ycenter for source
radius -- radius to measure the total number of counts on
"""
      exit(1)

   hdu = fits.open(sys.argv[1])
   xc = float(sys.argv[2])    
   yc = float(sys.argv[3])    
   radius = float(sys.argv[4])    

   #remove the background

   #determine the photmetry
   counts = aperture_phot(hdu[0].data, xc, yc, radius)

   #report the results
   print "No background Subtraction: ", counts

   
   

   #determine it with annulus background subtraction
   br1 =  radius+10 #these are arbitrary
   br2 =  radius+15
   counts1 = aperture_phot(hdu[0].data, xc, yc, br1)
   counts2 = aperture_phot(hdu[0].data, xc, yc, br2)
   #determine the background per pixels by normalizing the two 
   #annulli by their area
   bkgrd = (counts2 - counts1) / (np.pi * (br2**2 - br1**2))
   print "Annular background Subtraction: ",counts - bkgrd * np.pi * radius**2, bkgrd


   #determine the photometry with global background subtraction
   bkgrd = np.median(sigma_clip(hdu[0].data, sig=3, iters=5, cenfunc=np.median, varfunc=np.var)[0])
   gcounts = aperture_phot(hdu[0].data-bkgrd, xc, yc, radius)
   print "Global Background Subtract: ", gcounts, bkgrd
   
   #dispaly an aperture profile of the source and compare the average to the individual pixels
   r_arr = np.arange(1,radius+15)
   sb_arr = np.empty_like(r_arr)
   data = hdu[0].data-bkgrd
   for i,r in enumerate(r_arr): 
       sb_arr[i] = aperture_phot(data, xc, yc, r)


   plot(r_arr, sb_arr)
   show()

