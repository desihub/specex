#!/usr/bin/env python


import sys,string
import astropy.io.fits as pyfits

def usage() :
    print sys.argv[0], "-k key1 -k key2 ... file1 file2 ... (-hdu #)"
    sys.exit(0)


if len(sys.argv)<2 :
    usage()

filenames=[]
keys=[]
hdu=0

i=1
while i<len(sys.argv) :

    arg=sys.argv[i]
    if arg[0]!='-' :
        filenames.append(arg)
    else :
        if arg=="-k" :
            i+=1
            key=sys.argv[i]
            keys.append(key)

        elif arg=="-hdu" :
            i+=1
            hdu=string.atoi(sys.argv[i])
        
        else :
            print "unexpected arg",arg
            usage()
    i+=1

if len(filenames)==0 :
    print "no file given"
    usage()
    
if len(keys)==0 :
    print "no key, will dump header"

    for filename in filenames :
        print filename
        print "======================"
        h=pyfits.open(filename)
        header=h[hdu].header
        print header.tostring
        print ""
        h.close()
    sys.exit(0)
    
line="#"
for k in keys :
    line+=" "+k
line+=" filename"

print line
for filename in filenames :
    h=pyfits.open(filename)
    header=h[hdu].header
    line=""    
    for k in keys :
        try :
            val=header[k]
        except KeyError :
            val="None"
        line+=" "+str(val)
    line+=" "+filename
    print line
    h.close()

    
        
    
