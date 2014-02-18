#! /usr/bin/env python

"""
What:
Resize all the jpg images in a directory
All images will be placed inside a directory called "resized"
This directory must exist
Usage: ./resize.py
"""

import PIL.Image as PIL
import re, os, sys, urlparse

SIZE = 300,500
JPG = re.compile(".*\.(jpg|jpeg)", re.IGNORECASE)

def get_new_path(path):
    try:
        basedir = os.path.dirname(path) + '/resized/'
        if os.path.isdir(basedir):
            pass
        elif os.path.isfile(basedir):
            raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % basedir)
	else:
	    os.mkdir(basedir)
    except:
	print "Could not create the folder resized. Check privileges."
	exit(0)
    base, ext = os.path.splitext(os.path.basename(path))
    file = base + "_tn" + ext

    return urlparse.urljoin(basedir, file)

try:
    image_dir = sys.argv[1]
except:
    print "You must specify a directory"
    exit(0)

files = os.listdir(image_dir)
for file in files:
    if JPG.match(file):
        f = image_dir.rstrip("/") + "/" + file
        print f
        img = PIL.open(f)
        img.thumbnail(SIZE, PIL.ANTIALIAS)
        img.save(get_new_path(f))
