#!/usr/bin/python


import pyfits,sys

if len(sys.argv)<2 :
    print sys.argv[0],"table.fits"
    sys.exit(12);

filename=sys.argv[1]
# filename="psf-b1-00108382.fits"

hdulist=pyfits.open(filename)
hdulist.info()
hdu=hdulist[1]

print hdu.header
table = hdu.data
#print table
#print "COMMENT ------------------------------------------------------------------------"
print "--------------------------------------------------------------------------------"
cols = hdu.columns
print cols.names
print "--------------------------------------------------------------------------------"

nrows=table.shape[0]

for r in range(nrows) :
    line="%02d "%(r+1)
    line += "%s\t"%table.field("PARAM")[r]
    line += " [%d,%d]"%(table.field("WAVEMIN")[r],table.field("WAVEMAX")[r])
    
    coef = table.field("COEFF")[r]
    line += " %s"%str(coef.shape)
    line += " %f ..."%coef[0]
    print line


