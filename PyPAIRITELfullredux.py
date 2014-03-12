"""
PyPAIRITELfullredux.py reduces your PAIRITEL data. It extracts both 
long-read (exposure time 7.800 seconds) and short-read (exposure time 0.051 
seconds) data. It mosaics all the data and fits the mosaics with WCS header 
info.

Run the program in a directory with  your flats in a subdirectory named 
"./flats", and your bad pixel masks in a subdirectory named "./bpms". The 
program will copy everything to your local machine's /home/ptelreducer/scratch
and run the reductions
there. Once done, it will copy the reduction output directory to the place
specified with the -d argument.

The alignment software needs to be in a subdirectory named
"./alignment_software".

You will also need the following files:
    PyPAIRITELfullredux.py      -   This file, obviously
    pairitel_redux_j.sex        -   SExtractor configuration files for object
    pairitel_redux_h.sex            images production
    pairitel_redux_k.sex
    pairitel_redux.param        -   SExtractor output catalog parameters file
                                    for object images production
    gauss_3.0_7x7.conv          -   SExtractor filter for object image catalogs
    pairitel_cloud_detect.sex   -   SExtractor configuration file for cloudy
                                    triplestack image identification
    pairitel_cloud_detect.param -   SExtractor output catalog parameters file
                                    for cloudy triplestack image identification
    default.nnw                 -   SExtractor neural network file (necessary
                                    to determine object classificiation metric
                                    during cloudy image detection)
    anet.py                     -   Python code to use astrometry.net
    pairitel_redux.swarp        -   SWarp configuration file for mosacing

Dependencies:
    SExtractor
    SWarp
    WCSTools
    Python 2.5 or better
    scipy (python module)
    pyfits (python module)
    pyephem (python module)
    Ed Rosten's Astronomy Image Alignment Software (has it's own dependencies, 
        see below)

Modify the section USER-SPECIFIC PATHS AND VALUES below with the proper paths
to your SExtractor, SWarp, and sethead executables. Also, if you want to use a
python installation other than the default set by your environment variables,
you'll have to modify that, too.

The dependencies for Ed Rosten's code (./alignment_software) are:
    gawk (gnu awk, http://www.gnu.org/software/gawk/)
    TooN
    GVars3
    libCVD
    CVS (just to get them) - if you don't have CVS, check his site for 
        downloads: 
        http://mi.eng.cam.ac.uk/~er258/cvd/toon.html
        http://mi.eng.cam.ac.uk/~er258/cvd/index.html
        http://mi.eng.cam.ac.uk/~er258/cvd/gvars3.html
Download the latest versions from his CVS repository with the following 
commands.
TooN:
    cvs -z3 -d:pserver:anoncvs@cvs.savannah.nongnu.org:/cvsroot/toon co TooN
GVars3:
    cvs -z3 -d:pserver:anoncvs@cvs.savannah.nongnu.org:/cvsroot/libcvd co gvars3
libCVD:
    cvs -z3 -d:pserver:anoncvs@cvs.savannah.nongnu.org:/cvsroot/libcvd co libcvd
Install TooN first, then the others. To install these, just unzip, cd to the 
directory, and then do
    ./configure
    make
    sudo make install

You will need to edit your .bash_profile to contain the following lines:
    export DYLD_LIBRARY_PATH=/usr/local/lib/
    export LDFLAG=/usr/local/lib/

Once the dependencies for Ed's code are in place, cd alignment_software and run
    ./configure
    make
to compile the alignment software.

Once the program finishes it will dump the output into a subdirectory named
obs_string + "-reduction_output", where obs_string is the PAIRITEL database 
observation ID (ex: SN.175.3) of the first raw image in "./raw". This 
subdirectory contains sub-subdirectories for the mosaics, alignment files, and 
ancillary files. Depending on optional arguments it may also contain 
sub-subdirectories for the single-frame weightmaps, triplestacks, and reduced 
frames. All of these sub-subdirectories are labeled by PAIRITEL database object 
ID (ex: SN.175). The p0 frames are not used in the final mosaics because 
previous experience has shown them to exhibit abnormally high noise (either dark 
current or read noise). However, the program will dump the p0 triplestacks into 
the obj_string + "_triplestacks"/p0_triplestacks" subdirectory.

The triplestacks (save the p0 triplestacks) are aligned with Ed Rosten's
astronomical image alignment software. This results in pixel offsets wich are 
then applied to the triplestacks' WCS header info. So, the relative WCS info in 
the triplestacks is typically very good. The absolute WCS info in the 
triplestacks is generated with anet.py. It is also typically very 
good.

An explanation of the command-line arguments:
-r [path] REQUIRED specify path to directory containing the raw FITS files. scp 
is used to copy the files to your local /home/ptelreducer/scratch
directory. This means that the 
image files can be on your local machine or a distant machine (i.e., lyra) if 
you have passwordless ssh set up.
-o [obs_id] REQUIRED specify the PAIRITEL observation ID (e.g., PULSE.33.1) of 
the epoch you wish to reduce.
-d [path] REQUIRED specify the destination directory for the reduction products.
-s OPTIONAL use this to make final products of the short-read data.
-f OPTIONAL use this to make single-frame reduced and triplestacked final images
(the default is to delete these and only output final mosaics).
-p [fine_align_radius] OPTIONAL force alignment software to use telescope 
pointing offsets (culled from raw FITS headers) as initial guess and then search
within specified fine alignment radius to produce the final alignments.
-c OPTIONAL apply cloudy image rejection (default is to identify cloudy images
but not reject them from alignment or mosaicing).
-a [path to file] OPTIONAL use the specificed alignments file in lieu of running
alignment during reduction. This is most useful if you were unhappy with a 
previous reduction's alignment on the same raw dataset and modified the 
alignments file by hand. Also useful if you decide you wanted to use -f after
running the reduction already. This will skip the alignment process, which saves
a lot of computation time.

Debugging mode: If you set DEBUG = True in the USER-SPECIFIC PATHS AND VALUES
section then there will be an additional subdirectory created within the 
reduction_output directory called ObsID_debug. It will contain the midreduced
triplestacks, skarks, masks, and gaussian filtered object images. This can be
helpful if you want to figure out why a reduction failed or looks peculiar.

Below is an explanation of the actual reduction process.

The very first step is to subtract the 0.051 sec short-reads from the 7.851 sec
long-reads. This eliminates the bias in the long-reads (there is no way to 
remove it so precisely in the short-reads, instead we leave the bias in the 
short-read skarks and assume it doesn't change much over time). This step is 
done before making any midskarks or skarks and before running the images through
the reduction equation: reduced = (raw - skark)/flat.

The reduction proceeds by producing midskarks from the raw data. These are
not star-masked (the primary differences between these and the final skarks). 
Each rawfile has its own midskark (and skark) generated by the other time-local 
rawfiles. Further, each midskark (and skark) is particular to both band (J, H, 
or K) and position in the read cycle (0, 1, or 2). This second feature was added 
as an improvement over previous reduction codes which ignored read cycle 
position in skark generation. The default time window right now is +/- 5 minutes 
= 0.00347222 day. This can be altered in the code if necessary (it's in an "if" 
logic statement, and would need to be changefor both the code that makes 
midskarks and the code that makes skarks).

With midskarks in hand, the midreduction occurs. This is simply carrying out the
reduction equation as described above and transfering over FITS header values of
importance. Note that the exposure time and plate scale (SECPIX) are inserted
based on instrument specifications. The WCS information is inaccurate (although 
not horrible) and will be refit in the final mosaics. Also, the reduction
pipeline calculates and inserts the heliocentric julian date (of particular 
importance when studying variable sources over weeks or months).

The masks and dynamic_bpms are produced from the midreduced images. The masks
are applied during the average-combine to produce the final skarks so that the 
effects of any bright stars or galaxies are minimized in the skarks. If you see
"shadowing" around bright sources in the final mosaics it's a good bet that 
source flux made its way into the skarks.

The dynamic_bpms are generated by sigma-tagging the mid-reduced images
(-5, +3 sigma) so that only these very bright or very dim pixels are identified.
Then, a sum is made of these images to identify any pixels which are bright/dim
consistently regardless of dither position (thus, these pixels are actual 
defects on the detector, not bright stars). These sum images are turned into the 
band-read type specific dynamic_bpms and the union of these creates the 
dynamic_bpm which will be combined with the calibration bpms so that these 
pixels are de-weighted during the mosaicing.

The skarks are produced in nearly the same manner as the midskarks, except this 
time the masks are used to ignore flux from stars. Then, the final reduction 
proceeds in the same manner as the midreduction, except using the new, better 
skarks.

Triplestacks are produced by simply summing (pixel-wise) the three reads at each
dither position. Triplestacks of the midreduced data is used to create the 
alignment offsets (with Ed Rosten's Astronomy Image Alignment Software). The 
reasoning for using the midreduced triplestacks instead of the final 
triplestacks is that the final triplestacks have the potential to contain 
pixel regions with completely masked data. This is very rare and only occurs if
the masking was so intense and the dithers so poor that enough masks lined up to
create a void in the skark.

Once the triplestacks are aligned the "source image" triplestack is fit with
astrometry.net (or, at least we try). If this succeeds (common for denser star 
fields), then all the triplestacks are aligned well with absolute WCS. If this
fails, then the triplestacks have good relative alignment, but absolute WCS 
based on telescope pointing.

The mosaicing proceeds with SWarp and the final mosaics are fit with 
astrometry.net. It's quite rare that astrometry.net fails to fit a final mosaic.


Example Usage: 
    python PyPAIRITELfullredux.py 
     -r [user]@lyra.berkeley.edu:/Bloom/PAIRITEL-DATA/sem2009a/Dir2009-May-16/ 
     -o PULSE.33.1 
     -d /Users/CRK/Desktop/TCP_Repo/TCP/Software/Reduction_Pipeline_3

python2.5 PyPAIRITELfullredux.py -r /home/ptelreducer/PyPAIRITELfullredux_distribution/raw -o GRB.489.1 -d /home/ptelreducer/PyPAIRITELfullredux_distribution/ -f

python2.5 PyPAIRITELfullredux.py -r /home/ptelreducer/PyPAIRITELfullredux_distribution/raw -o PULSE.15.8 -d /home/ptelreducer/PyPAIRITELfullredux_distribution/ -f -s

python2.5 PyPAIRITELfullredux.py -r ptelreducer@lyra.berkeley.edu:/Bloom/PAIRITEL-DATA/sem2010a/Dir2010-Mar-26/ -o XRAYN.5.1 -d /home/ptelreducer/PyPAIRITELfullredux_distribution -f

On betsy computer:
python2.5 PyPAIRITELfullredux.py -r ptelreducer@lyra.berkeley.edu:/PAIRITEL/sem2010b/Dir2011-Jan-05/ -o PULSE.13.49 -d /home/ptelreducer/PyPAIRITELfullredux_distribution/reductions/PULSE/

python2.5 PyPAIRITELfullredux.py -r ptelreducer@lyra.berkeley.edu:/PAIRITEL/sem2010b/Dir2011-Jan-05/ -o PULSE.13.49 -d /home/ptelreducer/PyPAIRITELfullredux_distribution/reductions/PULSE/

python2.5 PyPAIRITELfullredux.py -r ptelreducer@lyra.berkeley.edu:/PAIRITEL/sem2010b/Dir2011-Jan-07/ -o GRB.567.9 -d /home/ptelreducer/PyPAIRITELfullredux_distribution/reductions/GRB/

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Author: Christopher Klein
Contact: cklein@astro.berkeley.edu
Date: began 2009/05/25
"""
#-------------------------------------------------------------------------------
# DECLARATION OF IMPORTS
import pyfits
from os import listdir, system, chdir, system, kill, popen
import sys
import ephem
from math import sin, cos, sqrt
from image import combine
from scipy import median
from scipy import sum
from scipy import std
from scipy import zeros
from scipy import where
from scipy import resize
from scipy import size
from scipy import clip
from scipy.ndimage import gaussian_filter
from processing import Pool
from processing import cpuCount
import operator
import signal
import threading
from optparse import OptionParser
#-------------------------------------------------------------------------------
# START TIMING
from time import time, sleep
# Begin program execution timing.
start_time = time()
#-------------------------------------------------------------------------------
# USER-SPECIFIC PATHS AND VALUES
# We'll be using SExtractor and SWarp. Define the proper bin executable here. 
sextractor_bin = "/home/ptelreducer/bin/sex"
swarp_bin = "/home/ptelreducer/bin/swarp"
# Same for SCAMP and MissFITS.
scamp_bin = "/home/ptelreducer/bin/scamp"
missfits_bin = "/home/ptelreducer/bin/missfits"
# We'll also be using some WCSTools.
sethead_bin = "/home/ptelreducer/bin/sethead"
# Tell the program which python installation you want to use.
python_bin = "/usr/bin/python2.5"
# For parallel processing we need the number of available CPUs. You could hard-
# code this value to something else, but it would result in non-optimal 
# performance.
numprocessors = cpuCount()
# Debugging mode outputs intermediate image files to help track down problems.
# It requires more time to run because of the additional hard disk reads/writes.
DEBUG = False
#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
# Quick function to divide up a large list into multiple small lists, 
# attempting to keep them all the same size. 
def split_seq(seq, size):
        newseq = []
        splitsize = 1.0/size*len(seq)
        for i in range(size):
            newseq.append(seq[int(round(i*splitsize)):
                int(round((i+1)*splitsize))])
        return newseq
# Simple scripted call to SExtractor to make object images.
def make_object_images(midreducedfilelist):
    for image_file in midreducedfilelist:
        if image_file[0] == "j":
            system(sextractor_bin + " midreduced/" + image_file + 
                " -c pairitel_redux_j.sex -WEIGHT_TYPE MAP_WEIGHT " + 
                "-WEIGHT_IMAGE bpms/jweightmap.fits " + 
                "-DETECT_MINAREA " + j_min_area + " " + 
                "-DETECT_THRESH " + j_thresh + " " +
                "-ANALYSIS_THRESH " + j_thresh + " " +
                "-CHECKIMAGE_NAME objectimages/objects-" + image_file)
        if image_file[0] == "h":
            system(sextractor_bin + " midreduced/" + image_file + 
                " -c pairitel_redux_h.sex -WEIGHT_TYPE MAP_WEIGHT " + 
                "-WEIGHT_IMAGE bpms/hweightmap.fits " + 
                "-DETECT_MINAREA " + h_min_area + " " + 
                "-DETECT_THRESH " + h_thresh + " " +
                "-ANALYSIS_THRESH " + h_thresh + " " +
                "-CHECKIMAGE_NAME objectimages/objects-" + image_file)
        # change pairitel_redux_k.sex to kill k band object detection and masking
        # pairitel_redux_k_nodetect.sex
        if image_file[0] == "k":
            system(sextractor_bin + " midreduced/" + image_file + 
                " -c pairitel_redux_k.sex -WEIGHT_TYPE MAP_WEIGHT " + 
                "-WEIGHT_IMAGE bpms/kweightmap.fits " + 
                "-DETECT_MINAREA " + k_min_area + " " + 
                "-DETECT_THRESH " + k_thresh + " " +
                "-ANALYSIS_THRESH " + k_thresh + " " +
                "-CHECKIMAGE_NAME objectimages/objects-" + image_file)
    return
# Simple function to write the command which will be used to call process2.
def make_process2_command(outputdir, source_imagefile, imagefile_list):
    process2_command = ("./process2 ../" + outputdir + " ../" + 
        source_imagefile)
    for imagefile in imagefile_list:
        process2_command = process2_command + " ../" + imagefile
    return process2_command
# Actually calling process2. The predefined function is necessary for 
# parallelization.
def run_process2(command):
    chdir("./alignment_software")
    system(command)
    chdir("..")
    return
# Simple function allowing parallelization of mosaicing with SWarp.
def run_swarp(command):
    system(command)
    return
# Simple function allowing parallelization of ncc_fine.
def run_ncc_fine(command):
    chdir("./alignment_software")
    for single_command in command:
        system(single_command)
    chdir("..")
    return
# Simple scripted call to SExtractor to make cloud_detect catalogs.
def detect_sources(triplestacks):
    for image_file in triplestacks:
        system(sextractor_bin + " " + image_file + 
            " -c pairitel_cloud_detect.sex -WEIGHT_TYPE MAP_WEIGHT " + 
            "-WEIGHT_IMAGE weightmaps/j_weightmap.fits " + 
            "-CATALOG_NAME cloud_cats/cloud_cat-" + image_file[13:-4] + "txt")
    return
# Median Absolute Deviation clipping for input list of numbers.
def mad_cliping(input_data):
    medval = median(input_data)
    sigma = 1.48 * median(abs(medval - input_data))
    high_sigma_clip_limit = medval + 3 * sigma
    low_sigma_clip_limit = medval - 3 * sigma
    clipped_data = input_data.clip(min=low_sigma_clip_limit, 
        max=high_sigma_clip_limit)
    new_medval = median(clipped_data)
    new_sigma = 1.48 * median(abs(medval - clipped_data))
    return clipped_data, new_medval, new_sigma
# All the rest is to get scamp to fit a single image (instead of anet.py)
def fit_wcs_scamp(image_filename, weightimage_filename):
    system(sextractor_bin + " " + image_filename + " " + 
        "-c pairitel_align.sex " + 
        "-WEIGHT_IMAGE " + weightimage_filename + " " +
        "-WEIGHT_TYPE MAP_WEIGHT " + 
        "-CATALOG_NAME " + image_filename.rstrip("fits") + "ldac")
    system(scamp_bin + " " + image_filename.rstrip("fits") + "ldac " + 
        "-c pairitel_align.scamp -ASTREF_CATALOG 2MASS")
    system(python_bin + " missfits.py " + image_filename)
    return
def test_wcs_scamp(image_filename):
    scamp_worked = False
    try:
        imagefile_hdulist = pyfits.open(image_filename)
        imagefile_header = imagefile_hdulist[0].header
    except:
        return scamp_worked
    try:
        scamp_comment = str(imagefile_header["PV1_1"])
        imagefile_hdulist.close()
    except:
        imagefile_hdulist.close()
        scamp_worked = False
        return scamp_worked
    try:
        if scamp_comment != "0.0":
            scamp_worked = True
    except:
        scamp_worked = False
    return scamp_worked
def kill_scamp(processname):
    sleep(5)
    for n in range(20):
        pid_found = False
        for line in popen("ps xa"):
            fields = line.split()
            if fields[4].split()[0] == processname:
                pid = fields[0]
                pid_found = True
                break
        if not pid_found:
            return
        else:
            sleep(1)
    try:
        kill(int(pid), signal.SIGHUP)
    except:
        return
    return
class scampThread(threading.Thread):
    def __init__ (self, image_filename, weightimage_filename):
        threading.Thread.__init__(self)
        self.image_filename = image_filename
        self.weightimage_filename = weightimage_filename
        self.status = -1
    def run(self):
        fit_wcs_scamp(self.image_filename, self.weightimage_filename)
class killThread(threading.Thread):
    def __init__ (self, scamp_bin):
        threading.Thread.__init__(self)
        self.scamp_bin = scamp_bin
    def run(self):
        global scamp_bin
        kill_scamp(self.scamp_bin)
def scamp_call(image_filename, weightimage_filename):
    threadlist = []
    current = scampThread(image_filename, weightimage_filename)
    threadlist.append(current)
    current.start()
    current = scampThread(image_filename.replace("j_long", "h_long"),
        weightimage_filename.replace("j_long", "h_long"))
    threadlist.append(current)
    current.start()
    current = scampThread(image_filename.replace("j_long", "k_long"),
        weightimage_filename.replace("j_long", "k_long"))
    threadlist.append(current)
    current.start()
    current = killThread(scamp_bin)
    threadlist.append(current)
    current.start()
    for job in threadlist:
        job.join()
    scamp_worked_j = False
    scamp_worked_j = test_wcs_scamp(image_filename)
    scamp_worked_h = False
    scamp_worked_h = test_wcs_scamp(image_filename.replace("j_long", "h_long"))
    scamp_worked_k = False
    scamp_worked_k = test_wcs_scamp(image_filename.replace("j_long", "k_long"))
    if scamp_worked_j and scamp_worked_h and scamp_worked_k:
        return True
    else:
        return False
#-------------------------------------------------------------------------------
# BEGIN MAIN PROGRAM
parser = OptionParser()

parser.add_option("-r", "--raw", "--raws", "--raw-directory", action="store",     
    type="string", dest="raw_directory", help=("path to directory containing " +
    "image[.fits] files"))
parser.add_option("-o", "--obs", "--obs-id", "--obs-string", action="store",     
    type="string", dest="obs_string", help=("observation ID string of epoch " + 
    "to be reduced"))
parser.add_option("-d", "--dest", "--dest-dir", action="store",     
    type="string", dest="destination_directory", 
    help="directory path to where the reduction output will be saved")   
parser.add_option("-s", "--short", "--shorts",
                  action="store_true", dest="do_short", default=False,
                  help="reduce short read data")
parser.add_option("-f", "--full", "--full-redux",
                  action="store_true", dest="do_full", default=False,
                  help=("full reduction to produce single images, " +              
                  "triplestacks, and mosaics (default is only final mosaics)"))
parser.add_option("-c", "--cloud", "--cloud-reject",
                  action="store_true", dest="do_cloud", default=False,
                  help=("do rejection of cloudy triplestacks before " + 
                  "alignment"))
parser.add_option("-p", "--pointing-alignment", action="store",     
    type="int", dest="pointing_radius", help=("use pointing offsets as coarse" + 
    " alignment and apply given radius as fine search mesh"))
parser.add_option("-a", "--align-file", action="store",     
    type="string", dest="alignment_file", 
    help="directory path of previously generated alignment file")

parser.add_option("--Ja", "--J-MinArea", action="store",     
    type="string", dest="j_min_area", default="120", 
    help="Minimum area in pixels for object detection by SExtractor during " + 
        "J band masking")
parser.add_option("--Ha", "--H-MinArea", action="store",     
    type="string", dest="h_min_area", default="100",
    help="Minimum area in pixels for object detection by SExtractor during " + 
        "H band masking")
parser.add_option("--Ka", "--K-MinArea", action="store",     
    type="string", dest="k_min_area", default="70",
    help="Minimum area in pixels for object detection by SExtractor during " + 
        "K band masking")

parser.add_option("--Jt", "--J-Thresh", action="store",     
    type="string", dest="j_thresh", default="4",
    help="Threshold (sigmas) for object detection/analysis by SExtractor during " + 
        "J band masking")
parser.add_option("--Ht", "--H-Thresh", action="store",     
    type="string", dest="h_thresh", default="4",
    help="Threshold (sigmas) for object detection/analysis by SExtractor during " + 
        "H band masking")
parser.add_option("--Kt", "--K-Thresh", action="store",     
    type="string", dest="k_thresh", default="4",
    help="Threshold (sigmas) for object detection/analysis by SExtractor during " + 
        "K band masking")

parser.add_option("--debug", action="store_true",     
    dest="DEBUG", default=False,
    help="Enable debugging mode")

(options, args) = parser.parse_args()
raw_directory = options.raw_directory
if not raw_directory:
    print "No raw directory specified, exiting."
    sys.exit()
obs_string = options.obs_string
if not obs_string:
    print "No obs_string specified, exiting."
    sys.exit()
dest_directory = options.destination_directory
if not dest_directory:
    print "No reduction destination directory specified, exiting."
    sys.exit()
do_short = options.do_short
do_full = options.do_full
do_cloud = options.do_cloud
pointing_radius = options.pointing_radius
if options.alignment_file == "None":
    alignment_file = False
else:
    alignment_file = options.alignment_file
j_min_area = options.j_min_area
h_min_area = options.h_min_area
k_min_area = options.k_min_area
j_thresh = options.j_thresh
h_thresh = options.h_thresh
k_thresh = options.k_thresh
DEBUG = options.DEBUG
system("rm -rf /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("mkdir /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("mkdir /home/ptelreducer/scratch/" + obs_string + "-reduction_directory/" + 
    obs_string + "_raw")
system("scp -P 10222 -q " + raw_directory + "/r*" + obs_string + "-*.fits /home/ptelreducer/scratch/" + 
    obs_string + "-reduction_directory/" + obs_string + "_raw")
system("scp -P 10222 -q " + raw_directory + "/Raw/r*" + obs_string + "-*.fits /home/ptelreducer/scratch/" + 
    obs_string + "-reduction_directory/" + obs_string + "_raw")
system("cp -r alignment_software /home/ptelreducer/scratch/" + obs_string + 
    "-reduction_directory/alignment_software")
system("cp -r flats /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp -r bpms /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp PyPAIRITELfullredux.py /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp pairitel_redux_j.sex /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp pairitel_redux_h.sex /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp pairitel_redux_k.sex /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp pairitel_redux_k_nodetect.sex /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp pairitel_redux.param /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp gauss_3.0_7x7.conv /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp pairitel_cloud_detect.sex /home/ptelreducer/scratch/" + obs_string + 
    "-reduction_directory")
system("cp pairitel_cloud_detect.param /home/ptelreducer/scratch/" + obs_string + 
    "-reduction_directory")
system("cp default.nnw /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp anet.py /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp anet_triplestack.py /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp autoastrometry.py /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp missfits.py /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp missfits2.py /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp pairitel_redux.swarp /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp pairitel_align.param /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp pairitel_align.scamp /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp pairitel_align.sex /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
system("cp gauss_2.0_3x3.conv /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
raw_directory = obs_string + "_raw"
full_raw_directory = "/home/ptelreducer/scratch/" + obs_string + "-reduction_directory/" + obs_string + "_raw"

# 11/23/2010 CRK: altered raw data clipping to only reduce the most recent
# series of observations on the given night. Observations can "stutter", so this
# will use only one run from the latest p0-0 onward. This is particularly 
# important with the new pointing check code that can cause an observation to 
# restart if it is found to have started at the wrong sky position. We really, 
# really don't want to reduce and coadd data mixed data (images from the wrong
# field with images from the right field). This will result in meaningless
# offsets from the alignment code and bad final mosaics.
rawfilelist = listdir(full_raw_directory)
rawfilelist2 = []
for rawfile in rawfilelist:
    if rawfile[-5:] == ".fits":
        rawfilelist2.append(rawfile)
obs_id = rawfilelist2[0].split("-")[4]

# Determine if we want to use the "new" version of rawfilelist creation or the old.
# The old version just lumps all images together and renumbers the p#-# to make 
# them sequential. It can handle "scheduling jitters". The new version ignores 
# all but the most recently obtained p0-0 and up set of images. If there was
# a "scheduling jitter" or an observation repeat, it only uses the most recent 
# observation. We will make the determination of which to use by measuring 
# the time between image acquistions. If any time separation is larger than 100 
# seconds, there was likely a pointing check forced restart. If all time 
# separations are less, there was probably just a scheduling jitter.

max_time_separation = 0
rawfilelist2.sort()
for n in range(len(rawfilelist2) - 1):
    year_a = rawfilelist2[n][1:5]
    month_a = rawfilelist2[n][6:9]
    day_a = rawfilelist2[n][10:12]
    hour_a = rawfilelist2[n][13:15]
    min_a = rawfilelist2[n][16:18]
    sec_a = rawfilelist2[n][19:21]
    year_b = rawfilelist2[n+1][1:5]
    month_b = rawfilelist2[n+1][6:9]
    day_b = rawfilelist2[n+1][10:12]
    hour_b = rawfilelist2[n+1][13:15]
    min_b = rawfilelist2[n+1][16:18]
    sec_b = rawfilelist2[n+1][19:21]
    if month_a == "Jan": month_a = "01"
    if month_a == "Feb": month_a = "02"
    if month_a == "Mar": month_a = "03"
    if month_a == "Apr": month_a = "04"
    if month_a == "May": month_a = "05"
    if month_a == "Jun": month_a = "06"
    if month_a == "Jul": month_a = "07"
    if month_a == "Aug": month_a = "08"
    if month_a == "Sep": month_a = "09"
    if month_a == "Oct": month_a = "10"
    if month_a == "Nov": month_a = "11"
    if month_a == "Dec": month_a = "12"
    if month_b == "Jan": month_b = "01"
    if month_b == "Feb": month_b = "02"
    if month_b == "Mar": month_b = "03"
    if month_b == "Apr": month_b = "04"
    if month_b == "May": month_b = "05"
    if month_b == "Jun": month_b = "06"
    if month_b == "Jul": month_b = "07"
    if month_b == "Aug": month_b = "08"
    if month_b == "Sep": month_b = "09"
    if month_b == "Oct": month_b = "10"
    if month_b == "Nov": month_b = "11"
    if month_b == "Dec": month_b = "12"
    first_time = ephem.date(year_a + "/" + month_a + "/" + day_a + " " + 
        hour_a + ":" + min_a + ":" + sec_a)
    second_time = ephem.date(year_b + "/" + month_b + "/" + day_b + " " + 
        hour_b + ":" + min_b + ":" + sec_b)
    separation_seconds = (second_time - first_time)*86400
    if separation_seconds > max_time_separation:
        max_time_separation = separation_seconds

    
use_new = True
if max_time_separation < 100: use_new = False

if use_new:
    rawfilelist2.sort(reverse=True)
    rawfilelist3 = []
    for item in rawfilelist2:
        rawfilelist3.append(item)
        if (item.split("-")[5] == "p0") and (item.split("-")[6] == "0.fits"):
            break
    rawfilelist3.sort()
    rawfilelist = rawfilelist3

if not use_new:
    rawfilelist2.sort()
    dither = 0
    read = 0
    for rawfile in rawfilelist2:
    #     first_part = rawfile.split("-p")[0] + "-p"
        first_part = (rawfile.split("-")[0] + "-" + rawfile.split("-")[1] + "-" + 
            rawfile.split("-")[2] + "-" + rawfile.split("-")[3] + "-" + obs_id + "-p")
        if read % 3 == 0:
            read = 0
        else:
            dither = dither - 1
        system("mv " + str(full_raw_directory) + "/" + rawfile + " " + 
            str(full_raw_directory) + "/" + first_part + str(dither) + "-" + 
            str(read) + ".fits")
        dither = dither + 1
        read = read + 1
    rawfilelist = listdir("/home/ptelreducer/scratch/" + obs_string + "-reduction_directory/" + raw_directory)

chdir("/home/ptelreducer/scratch/" + obs_string + "-reduction_directory")

# Create the error.log file. We'll append the PyFITS warning messages to this.
errorlogfile = open("./error.log", "w")
errorlogfile.write("PyPAIRITELfullredux.py error log.\n")
errorlogfile.close()
# Read in the raw images and build up the rawfiledictionary.
print "Now reading in raw images."
t1 = time()
# Create dictionary of raw data files. 
# These are files of the form r2009-May-16-10h38m18s-PULSE.74.1-p3-0.fits
# which are the raw output from the PAIRITEL 2MASS Survey cam.




rawfilelist_vetted = []
for filename in rawfilelist:
    timestring = filename.split("-")[3]
    time_seconds = (int(timestring.split("h")[0])*3600 + 
        int(timestring.split("h")[1].split("m")[0])*60 + 
        int(timestring.split("m")[1].split("s")[0]))
    rawfilelist_vetted.append([filename, timestring, time_seconds])
rawfilelist_vetted_sorted = sorted(rawfilelist_vetted, key=operator.itemgetter(0))
rawfilelist_vetted = rawfilelist_vetted_sorted
rawfilelist_vetted_2 = []
for n in range(len(rawfilelist_vetted)):
    if ((rawfilelist_vetted[n][0].split("-")[6].split(".fits")[0] == "0") and
        (rawfilelist_vetted[n+1][0].split("-")[6].split(".fits")[0] == "1") and
        (rawfilelist_vetted[n+2][0].split("-")[6].split(".fits")[0] == "2")):
        rawfilelist_vetted_2.append(rawfilelist_vetted[n])
    if ((rawfilelist_vetted[n][0].split("-")[6].split(".fits")[0] == "1") and
        (rawfilelist_vetted[n-1][0].split("-")[6].split(".fits")[0] == "0") and
        (rawfilelist_vetted[n+1][0].split("-")[6].split(".fits")[0] == "2")):
        rawfilelist_vetted_2.append(rawfilelist_vetted[n])
    if ((rawfilelist_vetted[n][0].split("-")[6].split(".fits")[0] == "2") and
        (rawfilelist_vetted[n-1][0].split("-")[6].split(".fits")[0] == "1") and
        (rawfilelist_vetted[n-2][0].split("-")[6].split(".fits")[0] == "0")):
        rawfilelist_vetted_2.append(rawfilelist_vetted[n])
rawfilelist_vetted_3 = []
rawfilelist_final = []
for n in range(len(rawfilelist_vetted_2)):
    skark_num = 0
    target_time = rawfilelist_vetted_2[n][2]
    for m in range(len(rawfilelist_vetted_2)):
        if (rawfilelist_vetted_2[n][0].split("-")[:6] == 
            rawfilelist_vetted_2[m][0].split("-")[:6]):
            continue
        elif ((abs(rawfilelist_vetted_2[m][2] - target_time) < 300) and 
            (rawfilelist_vetted_2[m][0].split("-")[6] == 
            rawfilelist_vetted_2[n][0].split("-")[6]) and 
            (rawfilelist_vetted_2[m][0].split("-")[5]) != "p0"):
            skark_num = skark_num + 1
    new_item = [rawfilelist_vetted_2[n][0], rawfilelist_vetted_2[n][1], 
        rawfilelist_vetted_2[n][2], skark_num]
    rawfilelist_vetted_3.append(new_item)
    if skark_num > 3:
        rawfilelist_final.append(rawfilelist_vetted_2[n][0])

rawfiledictionary = {}
for n in range(len(rawfilelist_final)):
    if rawfilelist_final[n][-5:] == ".fits":
        rawfiledictionary[rawfilelist_final[n]] = False
# Pluck the obs_string from the first raw file name. Ex: PULSE.74.3
obs_string = (rawfiledictionary.keys()[0].split("-")[4].split(".")[0] + "." + 
    rawfiledictionary.keys()[0].split("-")[4].split(".")[1] + "." + 
    rawfiledictionary.keys()[0].split("-")[4].split(".")[2])
# We now run through all the raw image files and fill in the rawfiledictionary.
# rawfiledictionary will have the file names as keys and another "subdictionary"
# as values. The subdictionaries will contain the importand header information
# and image data in the form of numarrays.
for rawfile in rawfiledictionary.keys():
    # Define (or clear out) the subdictionary.
    subdictionary = {}
    # Open the raw file as a Header Data Unit list. Pipe stdout again.
    stdout_orig = sys.stdout
    errorlogfile = open("./error.log", "a")
    sys.stdout = errorlogfile
    raw_hdulist = pyfits.open(raw_directory + "/" + rawfile)
    errorlogfile.close()
    sys.stdout = stdout_orig
    # Extract the fits header.
    rawheader = raw_hdulist[0].header
    # Pluck out DATE.
    observation_time = str(rawheader["DATE"])
    # Convert to ephem.date format and put in subdictionary.
    observation_time = (ephem.date(observation_time[:4] + "/" + 
        observation_time[5:7] + "/" + observation_time[8:10] + " " + 
        observation_time[11:13] + ":" + observation_time[14:16] + ":" + 
        observation_time[17:19]))
    subdictionary["observation_time"] = observation_time
    # Pluck out RA and DEC and put in subdictionary. These are in degrees.
    try:
        observation_ra = float(rawheader["RA"])
        subdictionary["observation_ra"] = observation_ra
    except(KeyError):
        observation_ra = float(rawheader["NEW_RA"])
        subdictionary["observation_ra"] = observation_ra
    try:
        observation_dec = float(rawheader["DEC"])
        subdictionary["observation_dec"] = observation_dec
    except(KeyError):
        observation_dec = float(rawheader["NEW_DEC"])
        subdictionary["observation_dec"] = observation_dec
    # Compute the observation time in Heliocentric Julian Date. First convert to 
    # julian date (jd) by adding the offset to ephem.date format.
    observation_time_jd = float(observation_time) + 2415020.0
    subdictionary["observation_time_jd"] = observation_time_jd
    # Convert observation_ra and observation_dec from degrees to float radians.
    observation_ra_radians = observation_ra * 0.01745329252
    observation_dec_radians = observation_dec * 0.01745329252
    # Calculate the Sun's position in the sky.
    sun = ephem.Sun()
    sun.compute(observation_time)
    sun_ra_radians = float(sun.a_ra)
    sun_dec_radians = float(sun.a_dec)
    # Calculate the Earth-Sun light travel time in days.
    earth_sun_light_travel_time = sun.earth_distance * 0.00577551833
    # Calculate the observation time in heliocentric julian date.
    observation_time_hjd = (observation_time_jd - earth_sun_light_travel_time * 
        (sin(observation_dec_radians) * sin(sun_ra_radians) + 
        cos(observation_dec_radians) * cos(sun_ra_radians) * 
        cos(observation_dec_radians - sun_dec_radians)))
    subdictionary["observation_time_hjd"] = observation_time_hjd
    # Pluck out a bunch of other important header information.
    subdictionary["pairitel_db_observation_id"] = rawheader["OBJECT"]
    subdictionary["target_name"] = rawheader["TRGTNAME"]
    subdictionary["strt_cpu"] = rawheader["STRT_CPU"]
    subdictionary["stop_cpu"] = rawheader["STOP_CPU"]
    try:
        subdictionary["principal_investigator"] = rawheader["PI"]
    except(KeyError):
        subdictionary["principal_investigator"] = "Unknown"
    try:
        subdictionary["airmass"] = rawheader["AIRMASS"]
    except(KeyError):
        subdictionary["airmass"] = "Not recorded"
    # Now extract the image data as numarrays. We immediately subtract the 
    # short-read from the long-read to eliminate the bias in the long-read.
    j_short_data = raw_hdulist[0].data[0]
    h_short_data = raw_hdulist[0].data[1]
    k_short_data = raw_hdulist[0].data[2]
    if do_short:
        subdictionary["j_short_data"] = j_short_data
        subdictionary["h_short_data"] = h_short_data
        subdictionary["k_short_data"] = k_short_data
    subdictionary["j_long_data"] = raw_hdulist[0].data[3] - j_short_data
    subdictionary["h_long_data"] = raw_hdulist[0].data[4] - h_short_data
    subdictionary["k_long_data"] = raw_hdulist[0].data[5] - k_short_data
    # Close up the raw fits file.
    raw_hdulist.close()
    # Now insert the subdictionary into rawfiledictionary with rawfie as key.
    rawfiledictionary[rawfile] = subdictionary
print "Finished reading in raw images, %f seconds required." % (time() - t1)
# Generate the midskarks and build up the midskarkdictionary.
print "Now making midskarks and performing midreduction."
t1 = time()
# Define the midreduceddictionary, which will be identical to rawfiledictionary
# in structure, but without header information (we don't need it).
midreduceddictionary = {}
# Read in the flats and create numarrays of their data.
jflat_data = pyfits.open("./flats/jflat.fits")[0].data
hflat_data = pyfits.open("./flats/hflat.fits")[0].data
kflat_data = pyfits.open("./flats/kflat.fits")[0].data
# This loop will compute the midskark data and build up midskarkdictionary.
for t_rawfile in rawfiledictionary.keys():
    # Clear out the subdictionary.
    subdictionary = {}
    # Store the target file's jd.
    target_jd = rawfiledictionary[t_rawfile]["observation_time_jd"]
    # Create (or clear out) lists to store midskark progenitor image data.
    if do_short:
        j_short_list = []
        h_short_list = []
        k_short_list = []
    j_long_list = []
    h_long_list = []
    k_long_list = []
    # Create (or clear out) the integer used to store the number of progenitor
    # images combined to produce each midskark.
    num_midskark_progenitors = 0
    # Then, we run through another loop of rawfiles to fill in the midskark 
    # progenitor image data lists.
    for p_rawfile in rawfiledictionary.keys():
        # Compare read positions ("0", "1", or "2"). We want the midskarks to 
        # be specific to read position.
        if t_rawfile.split("-")[6][0] == p_rawfile.split("-")[6][0]:
            # Store the potential progenitor file's jd.
            progenitor_jd = rawfiledictionary[p_rawfile]["observation_time_jd"]
            # Compare the target_jd and the progenitor_jd. Only progenitor 
            # images recorded within 5 minutes of the target image are used to 
            # make the midskark image. And, of course, we exclude the target 
            # file itself and the "p0" frames because their dark current is
            # inflated by the "down time" during the telescope slew.
            if (((target_jd - 0.00347222) < progenitor_jd) and
                ((target_jd + 0.00347222) > progenitor_jd) and 
                (target_jd != progenitor_jd) and
                (p_rawfile.split("-")[5] != "p0")):
                # Append the midskark progentor image data to the lists.
                if do_short:
                    j_short_list.append(
                        rawfiledictionary[p_rawfile]["j_short_data"])
                    h_short_list.append(
                        rawfiledictionary[p_rawfile]["h_short_data"])
                    k_short_list.append(
                        rawfiledictionary[p_rawfile]["k_short_data"])
                j_long_list.append(
                    rawfiledictionary[p_rawfile]["j_long_data"])
                h_long_list.append(
                    rawfiledictionary[p_rawfile]["h_long_data"])
                k_long_list.append(
                    rawfiledictionary[p_rawfile]["k_long_data"])
                # Increase the number of midskark progenitors for this target
                # image file.
                num_midskark_progenitors = num_midskark_progenitors + 1
    # Now make the midskark frames by average-combining.
    if do_short:
        j_short_ms = combine.average(j_short_list, nhigh=1, nlow=1)
        h_short_ms = combine.average(h_short_list, nhigh=1, nlow=1)
        k_short_ms = combine.average(k_short_list, nhigh=1, nlow=1)
    j_long_ms = combine.average(j_long_list, nhigh=1, nlow=1)
    h_long_ms = combine.average(h_long_list, nhigh=1, nlow=1)
    k_long_ms = combine.average(k_long_list, nhigh=1, nlow=1)
    # Now perform the reduction equation: ((raw - midskark) / flat) = midreduced
    if do_short:
        j_short_data = ((rawfiledictionary[t_rawfile]["j_short_data"] 
            - j_short_ms) / jflat_data)
        h_short_data = ((rawfiledictionary[t_rawfile]["h_short_data"] 
            - h_short_ms) / hflat_data)
        k_short_data = ((rawfiledictionary[t_rawfile]["k_short_data"] 
            - k_short_ms) / kflat_data)
    j_long_data = ((rawfiledictionary[t_rawfile]["j_long_data"] - j_long_ms)
        / jflat_data)
    h_long_data = ((rawfiledictionary[t_rawfile]["h_long_data"] - h_long_ms)
        / hflat_data)
    k_long_data = ((rawfiledictionary[t_rawfile]["k_long_data"] - k_long_ms)
        / kflat_data) 
    # Insert the midreduced data into the subdictionary.
    if do_short:
        subdictionary["j_short_data"] = j_short_data
        subdictionary["h_short_data"] = h_short_data
        subdictionary["k_short_data"] = k_short_data
    subdictionary["j_long_data"] = j_long_data
    subdictionary["h_long_data"] = h_long_data
    subdictionary["k_long_data"] = k_long_data
    midreduceddictionary["midreduced" + t_rawfile[1:]] = subdictionary
# We've got to write out the midreduced image files because SExtractor needs to
# read them in to create the object image files.
system("mkdir ./midreduced")
for midreducedfile in midreduceddictionary.keys():
    if do_short:
        output_hdu = pyfits.PrimaryHDU(
            midreduceddictionary[midreducedfile]["j_short_data"])
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdulist.writeto("./midreduced/j_short_" + midreducedfile)
        output_hdu = pyfits.PrimaryHDU(
            midreduceddictionary[midreducedfile]["h_short_data"])
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdulist.writeto("./midreduced/h_short_" + midreducedfile)
        output_hdu = pyfits.PrimaryHDU(
            midreduceddictionary[midreducedfile]["k_short_data"])
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdulist.writeto("./midreduced/k_short_" + midreducedfile)
    output_hdu = pyfits.PrimaryHDU(
        midreduceddictionary[midreducedfile]["j_long_data"])
    output_hdulist = pyfits.HDUList([output_hdu])
    output_hdulist.writeto("./midreduced/j_long_" + midreducedfile)
    output_hdu = pyfits.PrimaryHDU(
        midreduceddictionary[midreducedfile]["h_long_data"])
    output_hdulist = pyfits.HDUList([output_hdu])
    output_hdulist.writeto("./midreduced/h_long_" + midreducedfile)
    output_hdu = pyfits.PrimaryHDU(
        midreduceddictionary[midreducedfile]["k_long_data"])
    output_hdulist = pyfits.HDUList([output_hdu])
    output_hdulist.writeto("./midreduced/k_long_" + midreducedfile)
print ("Finished making midskarks and performing midreduction, " + 
    str(time() - t1) + " seconds required.")
# Use the midreduced image files and weightmaps (already written out) to feed 
# into SExtractor to generate object check images.
print "Now making object images."
t1 = time()
system("mkdir ./objectimages")
# Make a list of the midreduced file names.
midreducedfilelist = listdir("./midreduced")
midreducedfilelist2 = []
for n in range(len(midreducedfilelist)):
    if midreducedfilelist[n][-5:] == ".fits":
        midreducedfilelist2.append(midreducedfilelist[n])
midreducedfilelist = midreducedfilelist2
# Split the list for parallel processing.
split_midreducedfilelists = split_seq(midreducedfilelist, numprocessors)
# Run the object imaging making (through SExtractor) with parallel 
# processing.
p = Pool(numprocessors)
result = p.map_async(make_object_images, split_midreducedfilelists)
poolresult = result.get()
# Read in the object images and store them in a dictionary.
objectfilelist = listdir("./objectimages")
objectdictionary = {}
for n in range(len(objectfilelist)):
    if objectfilelist[n][-5:] == ".fits":
        objectdictionary[objectfilelist[n]] = False
for objectimage in objectdictionary.keys():
    objectdictionary[objectimage] = pyfits.open("./objectimages/" + 
        objectimage)[0].data
print ("Finished making object images, " + 
    str(time() - t1) + " seconds required.")
# Use the midreduced to make masks of saturated pixels and weightmaps.
print "Now making masks and weightmaps."
t1 = time()
# Create image_type_list - this will be used to loop through.
if do_short:
    image_type_list = ["j_short", "h_short", "k_short", 
        "j_long", "h_long", "k_long"]
if not do_short:
    image_type_list = ["j_long", "h_long", "k_long"]
dynamic_bpm_input_master_list = []
for midreducedfile in midreduceddictionary.keys():
    for image_type in image_type_list:
        midreduced_data = midreduceddictionary[midreducedfile][image_type + 
            "_data"]
        dynamic_bpm_input_master_list.append([midreduced_data, image_type, midreducedfile])


def make_dynamic_bpms(file_image_type_list):
    dynamic_bpm_flagged_data_dictionary = {}
    dynamic_bpm_flagged_data_dictionary_low = {}
    image_type_list = ["j_short", "h_short", "k_short", 
        "j_long", "h_long", "k_long"]
    for image_type in image_type_list:
        dynamic_bpm_flagged_data_dictionary[image_type] = zeros((256, 256))
        dynamic_bpm_flagged_data_dictionary_low[image_type] = zeros((256, 256))
    for entry in file_image_type_list:
        midreduced_data = entry[0]
        image_type = entry[1]
        midreducedfile = entry[2]
        midreduced_data_1d = resize(midreduced_data, size(midreduced_data))
        for n in range(3):
            midreduced_data_1d, medval, sigma = mad_cliping(midreduced_data_1d)
        # Define upper and lower limits for dynamically identified bad pixels.
        lower_limit = medval - 5 * sigma
        upper_limit = medval + 7 * sigma
        # Enforce upper and lower limits to create flagged image data.
        dynamic_bpm_flagged_data = (where(midreduced_data > upper_limit, 1, 0) + 
            where(midreduced_data < lower_limit, 1, 0))
        dynamic_bpm_flagged_data_low = where(midreduced_data < lower_limit, 1, 0)
        # Increment the dynamic_bpm_flagged_data.
        dynamic_bpm_flagged_data_dictionary[image_type] = (
            dynamic_bpm_flagged_data_dictionary[image_type] +             
            dynamic_bpm_flagged_data)
        dynamic_bpm_flagged_data_dictionary_low[image_type] = (
            dynamic_bpm_flagged_data_dictionary_low[image_type] +             
            dynamic_bpm_flagged_data_low)
    return [dynamic_bpm_flagged_data_dictionary, dynamic_bpm_flagged_data_dictionary_low]

# Split the list for parallel processing.
dynamic_bpm_input_master_list_split = split_seq(dynamic_bpm_input_master_list, numprocessors)

# Make masks using parallel processing.
p = Pool(numprocessors)
result = p.map_async(make_dynamic_bpms, dynamic_bpm_input_master_list_split)
poolresult = result.get()

# Create dynamic_bpm_flagged_data_dictionary and then fill it with zero-images.
dynamic_bpm_flagged_data_dictionary = {}
dynamic_bpm_flagged_data_dictionary_low = {}
for image_type in image_type_list:
    dynamic_bpm_flagged_data_dictionary[image_type] = zeros((256, 256))
    dynamic_bpm_flagged_data_dictionary_low[image_type] = zeros((256, 256))
    
for entry in poolresult:
    for image_type in image_type_list:
        dynamic_bpm_flagged_data_dictionary[image_type] = dynamic_bpm_flagged_data_dictionary[image_type] + entry[0][image_type]
        dynamic_bpm_flagged_data_dictionary_low[image_type] = dynamic_bpm_flagged_data_dictionary_low[image_type] + entry[1][image_type]
        
if DEBUG:
    system("rm -rf ./dynamic_bpm")
    system("mkdir ./dynamic_bpm")
    for image_type in image_type_list:
        output_hdu = pyfits.PrimaryHDU(
            dynamic_bpm_flagged_data_dictionary[image_type])
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdulist.writeto("./dynamic_bpm/" + str(image_type) + 
            "_high.fits")
        output_hdu = pyfits.PrimaryHDU(
            dynamic_bpm_flagged_data_dictionary_low[image_type])
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdulist.writeto("./dynamic_bpm/" + str(image_type) + 
            "_low.fits")
            
# We have to determine a value in the dynamic_bpm_flagged_data image below which
# a pixel is likely associated with a star and above which it is likely 
# associated with a dynamically identified broken pixel. The metric below is 
# pretty good. It says that any mask circles which overlap less than 2 out of 9 
# times are probably due to stars and permitted to stay. Anything else is weeded 
# out into the dynamic_bpm.
overlap_limit = len(midreduceddictionary.keys()) * (2.0/9.0)
# Enforce the overlap_limit.
for image_type in image_type_list:
    overlap_satisfied = where(dynamic_bpm_flagged_data_dictionary[image_type] > overlap_limit, 1, 0)
    low_satisfied = where(dynamic_bpm_flagged_data_dictionary_low[image_type] > overlap_limit, 1, 0)
    bpm_sum = overlap_satisfied + low_satisfied
    dynamic_bpm_flagged_data_dictionary[image_type] = where(bpm_sum > 0, 1, 0)
    
# Make weightmaps directory.
system("mkdir ./weightmaps")
# Read in and store the image data for each of the calibration bad pixel masks.
jbpm_hdulist = pyfits.open("./bpms/jbpm.fits")
jbpm_imagedata = jbpm_hdulist[0].data
jbpm_hdulist.close()
hbpm_hdulist = pyfits.open("./bpms/hbpm.fits")
hbpm_imagedata = hbpm_hdulist[0].data
hbpm_hdulist.close()
kbpm_hdulist = pyfits.open("./bpms/kbpm.fits")
kbpm_imagedata = kbpm_hdulist[0].data
kbpm_hdulist.close()
# Make the weightmap image data by combining the calibration bad pixel masks and
# the dynaimc_bpm_all_data.
jbpm_imagedata = where((jbpm_imagedata + where(( 
        dynamic_bpm_flagged_data_dictionary["j_long"]) > 0, 1, 0)) > 0, 0.0, 1.0)
hbpm_imagedata = where((hbpm_imagedata + where(( 
        dynamic_bpm_flagged_data_dictionary["h_long"]) > 0, 1, 0)) > 0, 0.0, 1.0)
kbpm_imagedata = where((kbpm_imagedata + where(( 
        dynamic_bpm_flagged_data_dictionary["k_long"]) > 0, 1, 0)) > 0, 0.0, 1.0)
# Write out the weightmaps.
output_hdu = pyfits.PrimaryHDU(jbpm_imagedata)
output_hdulist = pyfits.HDUList([output_hdu])
output_hdulist.writeto("./weightmaps/j_weightmap.fits")
output_hdu = pyfits.PrimaryHDU(hbpm_imagedata)
output_hdulist = pyfits.HDUList([output_hdu])
output_hdulist.writeto("./weightmaps/h_weightmap.fits")
output_hdu = pyfits.PrimaryHDU(kbpm_imagedata)
output_hdulist = pyfits.HDUList([output_hdu])
output_hdulist.writeto("./weightmaps/k_weightmap.fits")
# Create the maskdictionary.
maskdictionary = {}
# Then, fill in the maskdictionary with subdictionaries containing the mask 
# image data.
if DEBUG:
    system("rm -rf ./gaussianfiltered")
    system("mkdir ./gaussianfiltered")
    system("rm -rf ./sigmagaussian")
    system("mkdir ./sigmagaussian")
j_long_midreducedfilelist = []
for midreducedfile in midreducedfilelist:
    if midreducedfile[0] == "j" and midreducedfile[2] == "l":
            j_long_midreducedfilelist.append(midreducedfile)

midreduced_master_list = []
for n in range(len(j_long_midreducedfilelist)):
    midreducedfile = ("midreduced" + 
        j_long_midreducedfilelist[n].split("midreduced")[1])
    # image_type_list = ["j_short", "h_short", "k_short", "j_long", "h_long", "k_long"]
    # image_type_list = ["j_long", "h_long", "k_long"]
    for image_type in image_type_list:
        midreduced_data = (midreduceddictionary[midreducedfile][image_type + 
            "_data"])
        object_data = (objectdictionary["objects-" + image_type + "_" 
            + midreducedfile])
        midreduced_master_list.append([midreducedfile, image_type, 
            midreduced_data, object_data, DEBUG])

def make_masks(midreduced_master_list):
    mask_list = []
    for item in midreduced_master_list:
        midreducedfile = item[0]
        image_type = item[1]
        midreduced_data = item[2]
        object_data = item[3]
        DEBUG = item[4]
        mask_data = zeros((256, 256))
        if image_type == "j_long" or image_type == "j_short":
            midreduced_data_1d = resize(midreduced_data, size(midreduced_data))
            for n in range(3):
                midreduced_data_1d, medval, sigma = mad_cliping(
                    midreduced_data_1d)
            upper_limit = medval + 3.5 * sigma
            lower_limit = medval - 20.5 * sigma
            tagged_image0 = where(midreduced_data < lower_limit, 1, 0)
            tagged_image1 = where(midreduced_data > upper_limit, 1, 0)
            tagged_image2 = where(midreduced_data < -10000, 1, 0)
            sum_tagged = tagged_image0 + tagged_image1 + tagged_image2
            sigma_tagged_image = where((where(sum_tagged > 0, 1, 0) 
                - dynamic_bpm_flagged_data_dictionary[image_type]) > 0, 1, 0)
            blurred_sigma_tagged_image = zeros((256, 256))
            gaussian_filter(sigma_tagged_image, 0.5, 
                output=blurred_sigma_tagged_image)
            sigma_clipped_image = where(blurred_sigma_tagged_image > .01, 1, 0)
            if DEBUG:
                output_hdu = pyfits.PrimaryHDU(blurred_sigma_tagged_image)
                output_hdulist = pyfits.HDUList([output_hdu])
                output_hdulist.writeto("./sigmagaussian/sigmagaussian-" + 
                    image_type + "_" + midreducedfile)    
            object_tagged_image = where(object_data > 0, 20000/object_data, 0)
            gaussian_tagged_image = zeros((256, 256))
            gaussian_filter(object_tagged_image, 6, 
                output=gaussian_tagged_image)
            gaussian_tagged_image = 10*gaussian_tagged_image
            if DEBUG:
                output_hdu = pyfits.PrimaryHDU(gaussian_tagged_image)
                output_hdulist = pyfits.HDUList([output_hdu])
                output_hdulist.writeto("./gaussianfiltered/gaussianfiltered-" + 
                    image_type + "_" + midreducedfile)
            gaussian_clipped_image = where(gaussian_tagged_image > 5, 1, 0)
            mask_data = ((gaussian_clipped_image + sigma_clipped_image) > 0)
        if image_type == "h_long" or image_type == "h_short":
            midreduced_data_1d = resize(midreduced_data, size(midreduced_data))
            for n in range(3):
                midreduced_data_1d, medval, sigma = mad_cliping(
                    midreduced_data_1d)
            upper_limit = medval + 3.5 * sigma
            lower_limit = medval - 20.5 * sigma
            tagged_image0 = where(midreduced_data < lower_limit, 1, 0)
            tagged_image1 = where(midreduced_data > upper_limit, 1, 0)
            tagged_image2 = where(midreduced_data < -10000, 1, 0)
            sum_tagged = tagged_image0 + tagged_image1 + tagged_image2
            sigma_tagged_image = where((where(sum_tagged > 0, 1, 0) 
                - dynamic_bpm_flagged_data_dictionary[image_type]) > 0, 1, 0)
            blurred_sigma_tagged_image = zeros((256, 256))
            gaussian_filter(sigma_tagged_image, 0.5, 
                output=blurred_sigma_tagged_image)
            sigma_clipped_image = where(blurred_sigma_tagged_image > .01, 1, 0)
            if DEBUG:
                output_hdu = pyfits.PrimaryHDU(blurred_sigma_tagged_image)
                output_hdulist = pyfits.HDUList([output_hdu])
                output_hdulist.writeto("./sigmagaussian/sigmagaussian-" + 
                    image_type + "_" + midreducedfile)
            object_tagged_image = where(object_data > 0, 20000/object_data, 0)
            gaussian_tagged_image = zeros((256, 256))
            gaussian_filter(object_tagged_image, 6, 
                output=gaussian_tagged_image)
            gaussian_tagged_image = 10*gaussian_tagged_image
            if DEBUG:
                output_hdu = pyfits.PrimaryHDU(gaussian_tagged_image)
                output_hdulist = pyfits.HDUList([output_hdu])
                output_hdulist.writeto("./gaussianfiltered/gaussianfiltered-" + 
                    image_type + "_" + midreducedfile)
            gaussian_clipped_image = where(gaussian_tagged_image > 5, 1, 0)         
            mask_data = ((gaussian_clipped_image + sigma_clipped_image) > 0)
        if image_type == "k_long" or image_type == "k_short":
            midreduced_data_1d = resize(midreduced_data, size(midreduced_data))
            for n in range(3):
                midreduced_data_1d, medval, sigma = mad_cliping(
                    midreduced_data_1d)
            upper_limit = medval + 3.5 * sigma
            lower_limit = medval - 20.5 * sigma
            tagged_image0 = where(midreduced_data < lower_limit, 1, 0)
            tagged_image1 = where(midreduced_data > upper_limit, 1, 0)
            tagged_image2 = where(midreduced_data < -10000, 1, 0)
            sum_tagged = tagged_image0 + tagged_image1 + tagged_image2
            sigma_tagged_image = where((where(sum_tagged > 0, 1, 0) 
                - dynamic_bpm_flagged_data_dictionary[image_type]) > 0, 1, 0)
            blurred_sigma_tagged_image = zeros((256, 256))
            gaussian_filter(sigma_tagged_image, 0.5, 
                output=blurred_sigma_tagged_image)
            sigma_clipped_image = where(blurred_sigma_tagged_image > .01, 1, 0)
            if DEBUG:
                output_hdu = pyfits.PrimaryHDU(blurred_sigma_tagged_image)
                output_hdulist = pyfits.HDUList([output_hdu])
                output_hdulist.writeto("./sigmagaussian/sigmagaussian-" + 
                    image_type + "_" + midreducedfile)
            object_tagged_image = where(object_data > 0, 20000/object_data, 0)
            gaussian_tagged_image = zeros((256, 256))
            gaussian_filter(object_tagged_image, 6, 
                output=gaussian_tagged_image)
            gaussian_tagged_image = 10*gaussian_tagged_image
            if DEBUG:
                output_hdu = pyfits.PrimaryHDU(gaussian_tagged_image)
                output_hdulist = pyfits.HDUList([output_hdu])
                output_hdulist.writeto("./gaussianfiltered/gaussianfiltered-" + 
                    image_type + "_" + midreducedfile)
            gaussian_clipped_image = where(gaussian_tagged_image > 5, 1, 0)
            mask_data = ((gaussian_clipped_image + sigma_clipped_image) > 0)
        mask_list.append([midreducedfile, image_type, mask_data])
    return mask_list
                
# Split the list for parallel processing.
midreduced_master_list_split = split_seq(midreduced_master_list, numprocessors)

# Make masks using parallel processing.
p = Pool(numprocessors)
result = p.map_async(make_masks, midreduced_master_list_split)
poolresult = result.get()

mask_master_list = []
for m in range(len(poolresult)):
    for n in range(len(poolresult[m])):
        mask_master_list.append(poolresult[m][n])

mask_master_list_split = split_seq(mask_master_list, len(mask_master_list)/len(image_type_list))

for entry in mask_master_list_split:
    midreducedfile = entry[0][0]
    subdictionary = {}
    for n in range(len(image_type_list)):
        subdictionary[image_type_list[n] + "_data"] = entry[n][2]
    maskdictionary["mask" + midreducedfile[10:]] = subdictionary


if DEBUG:
    # Write out mask files.
    system("rm -rf ./masks")
    system("mkdir ./masks")
    for maskfile in maskdictionary.keys():
        if do_short:
            output_hdu = pyfits.PrimaryHDU(where(
                maskdictionary[maskfile]["j_short_data"], 1, 0))
            output_hdulist = pyfits.HDUList([output_hdu])
            output_hdulist.writeto("./masks/j_short_" + maskfile)
            output_hdu = pyfits.PrimaryHDU(where(
                maskdictionary[maskfile]["h_short_data"], 1, 0))
            output_hdulist = pyfits.HDUList([output_hdu])
            output_hdulist.writeto("./masks/h_short_" + maskfile)
            output_hdu = pyfits.PrimaryHDU(where(
                maskdictionary[maskfile]["k_short_data"], 1, 0))
            output_hdulist = pyfits.HDUList([output_hdu])
            output_hdulist.writeto("./masks/k_short_" + maskfile)
        output_hdu = pyfits.PrimaryHDU(where(
            maskdictionary[maskfile]["j_long_data"], 1, 0))
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdulist.writeto("./masks/j_long_" + maskfile)
        output_hdu = pyfits.PrimaryHDU(where(
            maskdictionary[maskfile]["h_long_data"], 1, 0))
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdulist.writeto("./masks/h_long_" + maskfile)
        output_hdu = pyfits.PrimaryHDU(where(
            maskdictionary[maskfile]["k_long_data"], 1, 0))
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdulist.writeto("./masks/k_long_" + maskfile)
print ("Finished making masks and weightmaps, " + 
    str(time() - t1) + " seconds required.")

print "Now making skarks and performing reduction."
if DEBUG:
    system("rm -rf ./skarks")
    system("mkdir ./skarks")
t1 = time()
# Define the reduceddictionary, which will be identical to rawfiledictionary in
# structure. Note that this time we do want to preserve header information, 
# since reduceddictionary will be used to write out the triplestacks.
reduceddictionary = {}
finalweightsdictionary = {}
# Recall that the flats have already been stored as numarrays.
# This loop will compute the skark data and build up skarkdictionary.
for t_rawfile in rawfiledictionary.keys():
    # Clear out the subdictionary.
    subdictionary = {}
    subdictionary2 = {}
    # Store the target file's jd.
    target_jd = rawfiledictionary[t_rawfile]["observation_time_jd"]
    # Create (or clear out) lists to store skark progenitor image data.
    if do_short:
        j_short_list = []
        h_short_list = []
        k_short_list = []
    j_long_list = []
    h_long_list = []
    k_long_list = []
    # Create (or clear out) lists to store the mask image data.
    if do_short:
        j_short_mask_list = []
        h_short_mask_list = []
        k_short_mask_list = []
    j_long_mask_list = []
    h_long_mask_list = []
    k_long_mask_list = []
    # Create (or clear out) the integer used to store the number of progenitor
    # images combined to produce each skark.
    num_skark_progenitors = 0
    # Then, we run through another loop of rawfiles to fill in the skark 
    # progenitor image data lists.
    for p_rawfile in rawfiledictionary.keys():
        # Compare read positions ("0", "1", or "2"). We want the skarks to 
        # be specific to read position.
        if t_rawfile.split("-")[6][0] == p_rawfile.split("-")[6][0]:
            # Store the potential progenitor file's jd.
            progenitor_jd = rawfiledictionary[p_rawfile]["observation_time_jd"]
            # Compare the target_jd and the progenitor_jd. Only progenitor 
            # images recorded within 5 minutes of the target image are used to 
            # make the skark image. And, of course, we exclude the target 
            # file itself and the "p0" frames because their dark current is
            # inflated by the "down time" during the telescope slew.
            if (((target_jd - 0.00347222) < progenitor_jd) and
                ((target_jd + 0.00347222) > progenitor_jd) and 
                (target_jd != progenitor_jd) and
                (p_rawfile.split("-")[5] != "p0")):
                # Compute and append the skark progentor image data to the 
                # lists. The objectdictionary makes subtraction of the object
                # flux relatively painless.
                # 090721: Commented object image subtraction from k band skarks.
                # This was causing a weird problem where the cores of stars 
                # in reduced images had pixel values < -50,000.
                if do_short:
                    j_short_list.append(
                        rawfiledictionary[p_rawfile]["j_short_data"])
                    h_short_list.append(
                        rawfiledictionary[p_rawfile]["h_short_data"])
                    k_short_list.append(
                        rawfiledictionary[p_rawfile]["k_short_data"])
                j_long_list.append(
                    rawfiledictionary[p_rawfile]["j_long_data"])
                h_long_list.append(
                    rawfiledictionary[p_rawfile]["h_long_data"])
                k_long_list.append(
                    rawfiledictionary[p_rawfile]["k_long_data"])
                # Collect the associated mask images into lists.
                if do_short:
                    j_short_mask_list.append(maskdictionary["mask" + 
                        p_rawfile[1:]]["j_short_data"])
                    h_short_mask_list.append(maskdictionary["mask" + 
                        p_rawfile[1:]]["h_short_data"])
                    k_short_mask_list.append(maskdictionary["mask" + 
                        p_rawfile[1:]]["k_short_data"])
                j_long_mask_list.append(maskdictionary["mask" + 
                    p_rawfile[1:]]["j_long_data"])
                h_long_mask_list.append(maskdictionary["mask" + 
                    p_rawfile[1:]]["h_long_data"])
                k_long_mask_list.append(maskdictionary["mask" + 
                    p_rawfile[1:]]["k_long_data"])
                # Increase the number of skark progenitors for this target
                # image file.
                num_skark_progenitors = num_skark_progenitors + 1
    # Now make the skark frames by average-combining. This time around we can
    # apply the masks we generated from the midreduced images.
    if do_short:
        j_short_s = combine.average(j_short_list, badmasks=j_short_mask_list,
            nhigh=1, nlow=1)
        h_short_s = combine.average(h_short_list, badmasks=h_short_mask_list, 
            nhigh=1, nlow=1)
        k_short_s = combine.average(k_short_list, badmasks=k_short_mask_list, 
            nhigh=1, nlow=1)
    j_long_s = combine.average(j_long_list, badmasks=j_long_mask_list, nhigh=1,
        nlow=1)
    h_long_s = combine.average(h_long_list, badmasks=h_long_mask_list, nhigh=1,
        nlow=1)
    k_long_s = combine.average(k_long_list, badmasks=k_long_mask_list, nhigh=1,
        nlow=1)
    num_skark_progenitors = num_skark_progenitors - 2
    if DEBUG:
        # Write out the skarks.
        if do_short:
            output_hdu = pyfits.PrimaryHDU(j_short_s)
            output_hdulist = pyfits.HDUList([output_hdu])
            output_hdr = output_hdulist[0].header
            output_hdr.update("NUMSKARK", num_skark_progenitors, 
                "Number of skark progenitor frames.")
            output_hdr.update("MEDSKARK", median(j_short_s), 
                "Median skark value.")
            output_hdulist.writeto("./skarks/j_short_skark" + t_rawfile[1:])
            output_hdu = pyfits.PrimaryHDU(h_short_s)
            output_hdulist = pyfits.HDUList([output_hdu])
            output_hdr = output_hdulist[0].header
            output_hdr.update("NUMSKARK", num_skark_progenitors, 
                "Number of skark progenitor frames.")
            output_hdr.update("MEDSKARK", median(h_short_s), 
                "Median skark value.")
            output_hdulist.writeto("./skarks/h_short_skark" + t_rawfile[1:])
            output_hdu = pyfits.PrimaryHDU(k_short_s)
            output_hdulist = pyfits.HDUList([output_hdu])
            output_hdr = output_hdulist[0].header
            output_hdr.update("NUMSKARK", num_skark_progenitors, 
                "Number of skark progenitor frames.")
            output_hdr.update("MEDSKARK", median(k_short_s), 
                "Median skark value.")
            output_hdulist.writeto("./skarks/k_short_skark" + t_rawfile[1:])
        output_hdu = pyfits.PrimaryHDU(j_long_s)
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdr = output_hdulist[0].header
        output_hdr.update("NUMSKARK", num_skark_progenitors, 
            "Number of skark progenitor frames.")
        output_hdr.update("MEDSKARK", median(j_long_s), 
            "Median skark value.")
        output_hdulist.writeto("./skarks/j_long_skark" + t_rawfile[1:])
        output_hdu = pyfits.PrimaryHDU(h_long_s)
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdr = output_hdulist[0].header
        output_hdr.update("NUMSKARK", num_skark_progenitors, 
            "Number of skark progenitor frames.")
        output_hdr.update("MEDSKARK", median(h_long_s), 
            "Median skark value.")
        output_hdulist.writeto("./skarks/h_long_skark" + t_rawfile[1:])
        output_hdu = pyfits.PrimaryHDU(k_long_s)
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdr = output_hdulist[0].header
        output_hdr.update("NUMSKARK", num_skark_progenitors, 
            "Number of skark progenitor frames.")
        output_hdr.update("MEDSKARK", median(k_long_s), 
            "Median skark value.")
        output_hdulist.writeto("./skarks/k_long_skark" + t_rawfile[1:])
    if do_short:
        skark_holes = where(j_short_s == 0, 0, 1)
        finalweights_data = where(skark_holes + jbpm_imagedata > 1, 1, 0)
        subdictionary2["j_short_data"] = finalweights_data
        skark_holes = where(h_short_s == 0, 0, 1)
        finalweights_data = where(skark_holes + hbpm_imagedata > 1, 1, 0)
        subdictionary2["h_short_data"] = finalweights_data
        skark_holes = where(k_short_s == 0, 0, 1)
        finalweights_data = where(skark_holes + kbpm_imagedata > 1, 1, 0)
        subdictionary2["k_short_data"] = finalweights_data
    skark_holes = where(j_long_s == 0, 0, 1)
    finalweights_data = where(skark_holes + jbpm_imagedata > 1, 1, 0)
    subdictionary2["j_long_data"] = finalweights_data
    skark_holes = where(h_long_s == 0, 0, 1)
    finalweights_data = where(skark_holes + hbpm_imagedata > 1, 1, 0)
    subdictionary2["h_long_data"] = finalweights_data
    skark_holes = where(k_long_s == 0, 0, 1)
    finalweights_data = where(skark_holes + kbpm_imagedata > 1, 1, 0)
    subdictionary2["k_long_data"] = finalweights_data
    finalweightsdictionary["finalweights" + t_rawfile[1:]] = subdictionary2
    # Now perform the reduction equation: ((raw - skark) / flat) = reduced
    if do_short:
        j_short_data = ((rawfiledictionary[t_rawfile]["j_short_data"] -
            j_short_s) / jflat_data)
        h_short_data = ((rawfiledictionary[t_rawfile]["h_short_data"] - 
            h_short_s) / hflat_data)
        k_short_data = ((rawfiledictionary[t_rawfile]["k_short_data"] - 
            k_short_s) / kflat_data)
    j_long_data = ((rawfiledictionary[t_rawfile]["j_long_data"] - j_long_s)
        / jflat_data)
    h_long_data = ((rawfiledictionary[t_rawfile]["h_long_data"] - h_long_s)
        / hflat_data)
    k_long_data = ((rawfiledictionary[t_rawfile]["k_long_data"] - k_long_s)
        / kflat_data) 
    # Insert the reduced data into the subdictionary.
    if do_short:
        subdictionary["j_short_data"] = j_short_data
        subdictionary["j_short_medskark"] = median(j_short_s)
        subdictionary["h_short_data"] = h_short_data
        subdictionary["h_short_medskark"] = median(h_short_s)
        subdictionary["k_short_data"] = k_short_data
        subdictionary["k_short_medskark"] = median(k_short_s)
    subdictionary["j_long_data"] = j_long_data
    subdictionary["j_long_medskark"] = median(j_long_s)
    subdictionary["h_long_data"] = h_long_data
    subdictionary["h_long_medskark"] = median(h_long_s)
    subdictionary["k_long_data"] = k_long_data
    subdictionary["k_long_medskark"] = median(k_long_s)
    reduceddictionary["reduced" + t_rawfile[1:]] = subdictionary
if DEBUG:
    system("rm -rf ./finalweights")
    system("mkdir ./finalweights")
    for finalweightsfile in finalweightsdictionary.keys():
        for image_type in image_type_list:
            output_hdu = pyfits.PrimaryHDU(
                finalweightsdictionary[finalweightsfile][image_type + "_data"])
            output_hdulist = pyfits.HDUList([output_hdu])
            output_hdulist.writeto("./finalweights/" + image_type + "_" + 
                finalweightsfile)
# Write out the reduced images. We can ignore WCS header info because we'll 
# copy it over from the triplestacks (after they are made and aligned).
system("mkdir ./reduced")
for reducedfile in reduceddictionary.keys():
    observation_time = str(rawfiledictionary[
        reducedfile.replace("educed", "")]["observation_time"])
    observation_time_hjd = rawfiledictionary[
        reducedfile.replace("educed", "")]["observation_time_hjd"]
    pairitel_db_observation_id = rawfiledictionary[
        reducedfile.replace("educed", "")]["pairitel_db_observation_id"]
    target_name = rawfiledictionary[
        reducedfile.replace("educed", "")]["target_name"]
    airmass = rawfiledictionary[
        reducedfile.replace("educed", "")]["airmass"]
    strt_cpu = rawfiledictionary[
        reducedfile.replace("educed", "")]["strt_cpu"]
    stop_cpu = rawfiledictionary[
        reducedfile.replace("educed", "")]["stop_cpu"]
    principal_investigator = rawfiledictionary[
        reducedfile.replace("educed", "")]["principal_investigator"]
    for image_type in image_type_list:
        output_hdu = pyfits.PrimaryHDU(
            reduceddictionary[reducedfile][image_type + "_data"])
        medskark_val = reduceddictionary[reducedfile][image_type + "_medskark"]
        output_hdulist = pyfits.HDUList([output_hdu])
        output_hdr = output_hdulist[0].header
        output_hdr.update("DATE-OBS", observation_time, "observation date, UTC")
        output_hdr.update("HJULDATE", observation_time_hjd, 
            "observation date, Heliocentric Julian Date")
        output_hdr.update("OBJECT",pairitel_db_observation_id , 
            "PAIRITEL observation ID")
        output_hdr.update("TRGTNAME",target_name , "target name")
        output_hdr.update("FILTER", image_type[0])
        if image_type[2:] == "long":
            output_hdr.update("EXPTIME", 7.800, 
                "exposure time, seconds")
        if image_type[2:] == "short":
            output_hdr.update("EXPTIME", 0.051, 
                "exposure time, seconds")
        output_hdr.update("GAIN", 10.5, "Nominal gain (e-/ADU)")
        output_hdr.update("AIRMASS", airmass, "secant z")
        output_hdr.update("MEDSKARK", medskark_val, "Median skark value.")
        output_hdr.update("STRT_CPU", strt_cpu)
        output_hdr.update("STOP_CPU", stop_cpu)
        output_hdr.update("PI", principal_investigator, 
            "principal investigator")
        output_hdr.update("INSTRUME", "2MASS Survey cam")
        output_hdr.update("OBSERVAT", "Mt. Hopkins")
        output_hdr.update("TELESCOP", "1.3m PAIRITEL")
        output_hdr.update("REDUX-SW", "PAIRITEL Pipeline v3.5", 
            "Version of reduction software")
        output_hdr.update("JDETAREA", j_min_area, 
            "Min detection area for J object masking")
        output_hdr.update("HDETAREA", h_min_area, 
            "Min detection area for H object masking")
        output_hdr.update("KDETAREA", k_min_area, 
            "Min detection area for K object masking")
        output_hdr.update("JTHRESH", j_thresh, 
            "Min detection threshold for J object masking")
        output_hdr.update("HTHRESH", h_thresh, 
            "Min detection threshold for H object masking")
        output_hdr.update("KTHRESH", k_thresh, 
            "Min detection threshold for K object masking")
        output_hdulist.writeto("./reduced/" + image_type + "_" + reducedfile)
print ("Finished making skarks and performing reduction, " + 
    str(time() - t1) + " seconds required.")
print "Now making midreduced triplestacks."
t1 = time()
system("mkdir ./triplestacks")
# Loop through reduceddictionary to median combine each dither position (p?-0,
# p?-1, and p?-2) and write out fits images with the header info we plucked out
# of the raw image files way back when we started.
for midreducedfile in midreduceddictionary.keys():
    # Store the filenamestring0, which is everything after "reduced".
    # Example: 2009-Jun-24-04h29m01s-SN.175.5-p4-0.fits
    filenamestring0 = midreducedfile.split("reduced")[1]
    if filenamestring0[-6] == "0":
        # Create filenamestring1 and filesnamestring2, which correspond to the 
        # other read positions.
        filenamestring1 = filenamestring0[:-6] + "1.fits"
        filenamestring2 = filenamestring0[:-6] + "2.fits"
        # Then, for each dither position combine the three read position frames
        # to create the triplestack data.
        for image_type in image_type_list:
#             triplestack_data = combine.median([
#                 midreduceddictionary["midreduced" + 
#                 filenamestring0][image_type + "_data"],
#                 midreduceddictionary["midreduced" + 
#                 filenamestring1][image_type + "_data"], 
#                 midreduceddictionary["midreduced" + 
#                 filenamestring2][image_type + "_data"]])
            triplestack_data = (
                midreduceddictionary["midreduced" + 
                filenamestring0][image_type + "_data"] + 
                midreduceddictionary["midreduced" + 
                filenamestring1][image_type + "_data"] + 
                midreduceddictionary["midreduced" + 
                filenamestring2][image_type + "_data"])
            # Write out the triplestack image. Transfer the header information
            # from the progenitor raw image files. Some of the header info is
            # known a priori for PAIRITEL, and thus didn't need to be plucked
            # from the raw files - we just write it in the triplestacks now.
            output_hdu = pyfits.PrimaryHDU(triplestack_data)
            output_hdulist = pyfits.HDUList([output_hdu])
            output_hdr = output_hdulist[0].header
            output_hdr.update("DATE-OBS", str(rawfiledictionary["r" + 
                filenamestring0]["observation_time"]), 
                "observation date, UTC")
            output_hdr.update("HJULDATE", rawfiledictionary["r" + 
                filenamestring0]["observation_time_hjd"], 
                "observation date, Heliocentric Julian Date")
            output_hdr.update("RA", rawfiledictionary["r" + 
                filenamestring0]["observation_ra"], 
                "telescope pointing, low accuracy")
            output_hdr.update("DEC", rawfiledictionary["r" + 
                filenamestring0]["observation_dec"], 
                "telescope pointing, low accuracy")
            output_hdr.update("EPOCH", 2000.0, 
                "Epoch")
            output_hdr.update("EQUINOX", 2000.0, 
                "Equinox")
            output_hdr.update("SECPIX", 2.0, 
                "Plate scale, arcseconds per pixel")
            output_hdr.update("CROTA2", 180.0, 
                "Image Twist +AXIS2 W of N, J2000 (deg)")
            output_hdr.update("OBJECT", rawfiledictionary["r" + 
                filenamestring0]["pairitel_db_observation_id"], 
                "PAIRITEL observation ID")
            output_hdr.update("TRGTNAME", rawfiledictionary["r" + 
                filenamestring0]["target_name"], 
                "target name")
            output_hdr.update("FILTER", image_type[0])
            if image_type[2:] == "long":
                output_hdr.update("EXPTIME", 23.400, 
                    "exposure time, seconds")
            if image_type[2:] == "short":
                output_hdr.update("EXPTIME", 0.153, 
                    "exposure time, seconds")
            output_hdr.update("GAIN", 10.5, "Nominal gain (e-/ADU)")
            output_hdr.update("AIRMASS", rawfiledictionary["r" + 
                filenamestring0]["airmass"], "secant z")
            output_hdr.update("STRT_CPU", rawfiledictionary["r" + 
                filenamestring0]["strt_cpu"])
            output_hdr.update("STOP_CPU", rawfiledictionary["r" + 
                filenamestring2]["stop_cpu"])
            output_hdr.update("STRT0000", rawfiledictionary["r" + 
                filenamestring0]["strt_cpu"])
            output_hdr.update("STOP0000", rawfiledictionary["r" + 
                filenamestring0]["stop_cpu"])
            output_hdr.update("STRT0001", rawfiledictionary["r" + 
                filenamestring1]["strt_cpu"])
            output_hdr.update("STOP0001", rawfiledictionary["r" + 
                filenamestring1]["stop_cpu"])
            output_hdr.update("STRT0002", rawfiledictionary["r" + 
                filenamestring2]["strt_cpu"])
            output_hdr.update("STOP0002", rawfiledictionary["r" + 
                filenamestring2]["stop_cpu"])
            output_hdr.update("PI", rawfiledictionary["r" + 
                filenamestring0]["principal_investigator"], 
                "principal investigator")
            output_hdr.update("INSTRUME", "2MASS Survey cam")
            output_hdr.update("OBSERVAT", "Mt. Hopkins")
            output_hdr.update("TELESCOP", "1.3m PAIRITEL")
            output_hdr.update("REDUX-SW", "PAIRITEL Pipeline v3.5", 
                "Version of reduction software")
            output_hdr.update("JDETAREA", j_min_area, 
                "Min detection area for J object masking")
            output_hdr.update("HDETAREA", h_min_area, 
                "Min detection area for H object masking")
            output_hdr.update("KDETAREA", k_min_area, 
                "Min detection area for K object masking")
            output_hdr.update("JTHRESH", j_thresh, 
                "Min detection threshold for J object masking")
            output_hdr.update("HTHRESH", h_thresh, 
                "Min detection threshold for H object masking")
            output_hdr.update("KTHRESH", k_thresh, 
                "Min detection threshold for K object masking")
            # Finally, just write out the image file.
            output_hdulist.writeto("./triplestacks/" + image_type + 
                "_triplestack" + filenamestring0[:-7] + ".fits")
print "Finished midreduced triplestacks, %f seconds required." % (time() - t1)
# Run image alignment.
print "Now detecting cloudy images and aligning the clear images."
t1 = time()
# This is the directory that will store SExtractor output catalogs which will 
# then be used to determine which images are cloudy.
system("mkdir cloud_cats")
# Make a list of all j_long triplestacks.
triplestackfilelist = listdir("./triplestacks")
j_long_filelist = []
for n in range(len(triplestackfilelist)):
    if triplestackfilelist[n][-5:] == ".fits":
        if triplestackfilelist[n][:6] == "j_long":
            j_long_filelist.append("triplestacks/" + triplestackfilelist[n])
# Split the j_long triplestack list into mutliple sublists for parallel 
# processing.
split_filelists = split_seq(j_long_filelist, numprocessors)
# Run detect_sources in parallel (calls SExtractor to make the cloud cats).
p = Pool(numprocessors)
result = p.map_async(detect_sources, split_filelists)
poolresult = result.get()
# Make a list of the cloud cats.
cloud_cat_list = []
for image_file in (j_long_filelist):
    cloud_cat_list.append("cloud_cats/cloud_cat-" + image_file[13:-4] + "txt")
# We analyze the cloud cats to count the number of good sources (detections with
# no abnormal SExtractor flags and > 0.68 star classification metric). We also
# keep track of the number of bad sources.
good_sources_list = []
good_sources_num_list = []
for catalog in cloud_cat_list:
    current_file = file(catalog, "r")
    num_good_sources = 0
    num_bad_sources = 0
    for line in current_file:
        if (int(line.split()[0]) == 0) and (float(line.split()[1]) > 0.68):
            num_good_sources = num_good_sources + 1
        else:
            num_bad_sources = num_bad_sources + 1
    current_file.close()
    good_sources_num_list.append(num_good_sources)
    good_sources_list.append(["triplestacks/" + catalog[21:-3] + "fits", 
        num_good_sources, num_bad_sources])
# Sort the list of images with good sources by the number of good sources.
good_sources_list_sorted = sorted(good_sources_list, key=operator.itemgetter(1))
# Calculate the median number of good sources and its standard deviation.
median_good_sources = median(good_sources_num_list)
stdev_good_sources = std(good_sources_num_list)
# Finally, apply these to determine which images are "cloudy" and which are 
# "clear".
all_triplestacks_list = []
cloudy_triplestacks_list = []
clear_triplestacks_list = []
for item in good_sources_list_sorted:
    # "Cloudy" is defined as having less than 3 sigma fewer good sources than
    # the median or having less than 2 good sources.
    if ((item[1] < (median_good_sources - 3 * stdev_good_sources)) or 
        (item[1] < 2)):
        item.append("cloudy")
        cloudy_triplestacks_list.append(item)
        all_triplestacks_list.append(item)
    else:
        item.append("clear")
        clear_triplestacks_list.append(item)
        all_triplestacks_list.append(item)
cloudy_clear_file = file("cloud_rejection.txt", "w")
for item in all_triplestacks_list:
    cloudy_clear_file.write(item[0] + " " + str(item[1]) + " " + str(item[2]) +
        " " + item[3] + "\n")
cloudy_clear_file.close()
# Sort the clear triplestacks list by number of good sources.
clear_triplestacks_list_sorted = sorted(clear_triplestacks_list, 
    key=operator.itemgetter(1))
cloudy_triplestacks_list_sorted = sorted(cloudy_triplestacks_list, 
    key=operator.itemgetter(1))
# Clear out the j_long_filelist; we're about to rewrite it with the clear 
# images.
j_long_filelist = []
# Define the source reference image (the image off of which the alignments will
# be calculated). Try to pick an image with lots of sources. But, to prevent
# possible problems due to possible reduction errors (artificial sources 
# produced by too many overlaping badmasks in the skarks), we don't select the
# image with the absolute greatest number of good sources (instead, we go with
# one with slightly fewer). Also, it's reasonable to make sure that the p0 frame
# is not used as the source image, since it won't be stacked.
try:
    if clear_triplestacks_list_sorted[-1][0].split("-p")[1][0] != "0":
        source_imagefile = clear_triplestacks_list_sorted[-1][0]
    else:
        source_imagefile = clear_triplestacks_list_sorted[-2][0]
except(IndexError):
    print "ERROR: No triplestacks passed the cloudy image filter, aborting."
    if cloudy_triplestacks_list_sorted[-1][0].split("-p")[1][0] != "0":
        source_imagefile = cloudy_triplestacks_list_sorted[-1][0]
    else:
        source_imagefile = cloudy_triplestacks_list_sorted[-2][0]
#     sys.exit()
# Now we can fill in the j_long_filelist, making sure to leave out the 
# source_imagefile because we don't need to calculate offsets for it.
for n in range(len(clear_triplestacks_list_sorted)):
    if clear_triplestacks_list_sorted[n][0] != source_imagefile:
        j_long_filelist.append(clear_triplestacks_list_sorted[n][0])
if not do_cloud:
    for n in range(len(cloudy_triplestacks_list_sorted)):
        if cloudy_triplestacks_list_sorted[n][0] != source_imagefile:
            j_long_filelist.append(cloudy_triplestacks_list_sorted[n][0])
j_long_filelist.sort()
# The final step in the cloudy image identification process is to move all the 
# cloudy triplestacks and reduced images into subdirectories.
system("mkdir triplestacks/cloudy")
system("mkdir reduced/cloudy")
if do_cloud:
    for cloudy_triplestack in cloudy_triplestacks_list:
        image_string = (cloudy_triplestack[0].split("j_long_")[0] + "*" + 
            cloudy_triplestack[0].split("j_long_")[1])
        system("mv " + image_string + " triplestacks/cloudy")
        image_string = image_string.replace("triplestacks", "reduced")
        image_string = image_string.replace("triplestack", "reduced")
        system("mv " + image_string[:-5] + "* reduced/cloudy")
    if len(cloudy_triplestacks_list) != 0:
        print "Detected %i cloudy triplestacks." % len(cloudy_triplestacks_list)
        print "Moving cloudy triplestacks and reduced images to subdirectory"
    else:
        print "Detected no cloudy triplestacks."
if not do_cloud:
    if len(cloudy_triplestacks_list) != 0:
        print "Detected %i cloudy triplestacks." % len(cloudy_triplestacks_list)
        print ("Cloud rejection is disabled by default, so not excluding " + 
            "cloudy triplestacks.")
    else:
        print "Detected no cloudy triplestacks."
# Pull RA and DEC from raw files to create estimated alignment offsets.
rawsourcefile = "r" + source_imagefile.strip("triplestacks/j_long_triplestack")
rawsourcefile = rawsourcefile[:-2] + "-1.fits"
rawsource_ra = rawfiledictionary[rawsourcefile]["observation_ra"]
rawsource_dec = rawfiledictionary[rawsourcefile]["observation_dec"]
ra_adjust_factor = cos(rawsource_dec * 0.01745329252)
pointingalignmentsfile = file("pointing_alignments.txt", "w")
j_long_full_filelist = [source_imagefile]
for filename in j_long_filelist:
    j_long_full_filelist.append(filename)
j_long_full_filelist.sort()
pointing_alignments_list = []
for j_long_triplestack in j_long_full_filelist:
    rawfile = "r" + j_long_triplestack.strip("triplestacks/j_long_triplestack")
    rawfile = rawfile[:-2] + "-1.fits"
    x_offset_deg = ((rawfiledictionary[rawfile]["observation_ra"] 
        - rawsource_ra) * ra_adjust_factor)
    x_offset_pix = x_offset_deg * 1800
    y_offset_deg = rawfiledictionary[rawfile]["observation_dec"] - rawsource_dec
    y_offset_pix = y_offset_deg * 1800
    pointingalignmentsfile.write(str(x_offset_pix) + " " + str(y_offset_pix) + 
        " " + j_long_triplestack + "\n")
    pointing_alignments_list.append([x_offset_pix, y_offset_pix, 
        j_long_triplestack])
pointingalignmentsfile.close()
if not alignment_file and not pointing_radius:
    # Now we can move on to calculating the alignments of the clear 
    # triplestacks. Split the triplestack list into mutliple sublists for 
    # parallel processing.
    print"Now running alignment software."
    split_filelists = split_seq(j_long_filelist, numprocessors)
    # Create the process2 commands.
    process2_command_list = []
    for num in range(numprocessors):
        process2_command_list.append(make_process2_command("process2output_"
            + str(num), source_imagefile, split_filelists[num]))
    # Run the process2 commands in parallel.
    p = Pool(numprocessors)
    result = p.map_async(run_process2, process2_command_list)
    poolresult = result.get()
    # Create the final_alignments_file for j band images, which will be made 
    # from the combination of the partial alignments files created by process2.
    final_alignments_list = []
    for num in range(numprocessors):
        try:
            partial_alignments_file = file("process2output_" + str(num) + 
                "/alignments.txt",  "r")
            if num == 0:
                for line in partial_alignments_file:
                    modified_line = line.replace("../", "").rstrip() + "\n"
                    final_alignments_list.append([modified_line.split()[0], 
                        modified_line.split()[1], modified_line.split()[2]])
            else:
                for line in partial_alignments_file:
                    if line.rstrip().split("../")[1] != source_imagefile:
                        modified_line = line.replace("../", "").rstrip() + "\n"
                        final_alignments_list.append([modified_line.split()[0], 
                            modified_line.split()[1], modified_line.split()[2]])
        except:
            continue
        partial_alignments_file.close()
    final_alignments_list_sorted = sorted(final_alignments_list, 
        key=operator.itemgetter(2))
    final_alignments_file = file("final_alignments_j.txt", "w")
    for item in final_alignments_list_sorted:
        final_alignments_file.write(item[0] + " " + item[1] + " " + 
            item[2] + "\n")
    final_alignments_file.close()
    
    
    
    
if not alignment_file and pointing_radius:
    mosaic_parameters_file = file("alignment_software/mosaic.cfg", "r")
    mosaic_parameters_file2 = file("alignment_software/mosaic2.cfg", "w")
    for line in mosaic_parameters_file:
        if line.split("=")[0] == "fine.mask.search_radius":
            line = "fine.mask.search_radius=" + str(pointing_radius) + "\n"
        mosaic_parameters_file2.write(line)
    mosaic_parameters_file.close()
    mosaic_parameters_file2.close()
    system("mv alignment_software/mosaic2.cfg alignment_software/mosaic.cfg")
        
    coarse_alignments_file = file("pointing_alignments.txt", "r")
    ncc_fine_commands = []
    for line in coarse_alignments_file:
        if line.split()[2] == source_imagefile:
            continue
        else:
            ncc_fine_commands.append('./ncc_fine --source ../' + 
                source_imagefile + ' --image ../' + 
                line.split()[2] + 
                ' --initial_offset "' + line.split()[0] + ' ' + 
                line.split()[1] + 
                ' " > ../j_alignments/fine/' + source_imagefile.split('/')[1] + 
                ',' + line.split()[2].split('/')[1] + 
                ',txt; cut -d\  -f 2,3 ../j_alignments/fine/' + 
                source_imagefile.split('/')[1] +',' + 
                line.split()[2].split('/')[1] + 
                ',txt | ./weighted_mean --start "' + line.split()[0] + ' ' + 
                line.split()[1] + ' " > ../j_alignments/final/' + 
                source_imagefile.split('/')[1] +',' + 
                line.split()[2].split('/')[1] + 
                ',txt; echo ' + line.split('-p')[1][:-6])
    coarse_alignments_file.close()        
    coarse_alignments_file = file("pointing_alignments.txt", "r")
    system("mkdir j_alignments")
    system("mkdir j_alignments/fine")
    system("mkdir j_alignments/final")
    split_ncc_fine_commands = split_seq(ncc_fine_commands, numprocessors)
    p = Pool(numprocessors)
    result = p.map_async(run_ncc_fine, split_ncc_fine_commands)
    poolresult = result.get()
    final_alignments_list = []
    alignment_file_j = file("final_alignments_j.txt", "w")
    for line in coarse_alignments_file:
        if line.split()[2] == source_imagefile:
            alignment_file_j.write("0 0 " + source_imagefile + "\n")
        else:
            j_final_alignment = file("j_alignments/final/" + 
                source_imagefile.split("/")[1] + "," + 
                line.split()[2].split("/")[1] + 
                ",txt", "r")
            for line2 in j_final_alignment:
                j_offsets = line2.rstrip()
            j_final_alignment.close()
            alignment_file_j.write(j_offsets + " " + 
                line.split()[2] + "\n")
    alignment_file_j.close()
    coarse_alignments_file.close()
    alignment_file_j = file("final_alignments_j.txt", "r")
    for line in alignment_file_j:
        final_alignments_list.append([line.split()[0], 
            line.split()[1], line.split()[2]])
    alignment_file_j.close()
    final_alignments_list_sorted = sorted(final_alignments_list, 
        key=operator.itemgetter(2))
    final_alignments_file = file("final_alignments_j.txt", "w")
    for item in final_alignments_list_sorted:
        final_alignments_file.write(item[0] + " " + item[1] + " " + 
            item[2] + "\n")
    final_alignments_file.close()



if alignment_file:
    system("cp " + alignment_file + " ./final_alignments_j.txt")
    final_alignments_list_sorted = []
    final_alignments_file = file("final_alignments_j.txt", "r")
    for line in final_alignments_file:
        final_alignments_list_sorted.append(line)
print "Finished alignments, %f seconds required." % (time() - t1)
print "Now making triplestacks."
t1 = time()
if DEBUG:
    system("mv ./triplestackweights ./midreduced_triplestackweights")
    system("mkdir ./triplestackweights")
    system("mv ./triplestacks ./midreduced_triplestacks")
    system("mkdir ./triplestacks")

else:    
    system("rm -rf ./triplestackweights")
    system("mkdir ./triplestackweights")
    system("rm -rf ./triplestacks")
    system("mkdir ./triplestacks")
reducedfilelist_clear = listdir("./reduced")
reducedfilelist_cloudy = listdir("./reduced/cloudy")
j_long_reducedfilelist = []
for reducedfile in reducedfilelist_clear:
    if (reducedfile[0] == "j" and reducedfile[2] == "l" and 
        reducedfile[-5:] == ".fits"):
        j_long_reducedfilelist.append(reducedfile)
for reducedfile in reducedfilelist_cloudy:
    if (reducedfile[0] == "j" and reducedfile[2] == "l" and 
        reducedfile[-5:] == ".fits"):
        j_long_reducedfilelist.append(reducedfile)
# Loop through reduceddictionary to median combine each dither position (p?-0,
# p?-1, and p?-2) and write out fits images with the header info we plucked out
# of the raw image files way back when we started.
for n in range(len(j_long_reducedfilelist)):
    if j_long_reducedfilelist[n][-6] == "0":
        reducedfile = j_long_reducedfilelist[n]
        # Store the filenamestring0, which is everything after "reduced".
        # Example: 2009-Jun-24-04h29m01s-SN.175.5-p4-0.fits
        filenamestring0 = reducedfile.split("reduced")[1]
        # Create filenamestring1 and filesnamestring2, which correspond to the 
        # other read positions.
        filenamestring1 = filenamestring0[:-6] + "1.fits"
        filenamestring2 = filenamestring0[:-6] + "2.fits"
        # Then, for each dither position combine the three read position frames
        # to create the triplestack data.
        for image_type in image_type_list:
#             triplestack_data = combine.median([
#                 reduceddictionary["reduced" + 
#                 filenamestring0][image_type + "_data"],
#                 reduceddictionary["reduced" + 
#                 filenamestring1][image_type + "_data"], 
#                 reduceddictionary["reduced" + 
#                 filenamestring2][image_type + "_data"]])
            triplestack_data = (
                reduceddictionary["reduced" + 
                filenamestring0][image_type + "_data"] + 
                reduceddictionary["reduced" + 
                filenamestring1][image_type + "_data"] + 
                reduceddictionary["reduced" + 
                filenamestring2][image_type + "_data"])
            weight_data = (where(
                finalweightsdictionary["finalweights" + 
                filenamestring0][image_type + "_data"] + 
                finalweightsdictionary["finalweights" + 
                filenamestring1][image_type + "_data"] + 
                finalweightsdictionary["finalweights" + 
                filenamestring2][image_type + "_data"] == 3, 1.0, 0.0))
            medskark_data = (
                reduceddictionary["reduced" + 
                filenamestring0][image_type + "_medskark"] + 
                reduceddictionary["reduced" + 
                filenamestring1][image_type + "_medskark"] + 
                reduceddictionary["reduced" + 
                filenamestring2][image_type + "_medskark"])
            output_hdu = pyfits.PrimaryHDU(weight_data)
            output_hdulist = pyfits.HDUList([output_hdu])
            output_hdr = output_hdulist[0].header
            output_hdulist.writeto("./triplestackweights/" + image_type + 
                "_triplestackweightmap" + filenamestring0[:-7] + ".fits")
            triplestack_data = triplestack_data * weight_data
            # Write out the triplestack image. Transfer the header information
            # from the progenitor raw image files. Some of the header info is
            # known a priori for PAIRITEL, and thus didn't need to be plucked
            # from the raw files - we just write it in the triplestacks now.
            output_hdu = pyfits.PrimaryHDU(triplestack_data)
            output_hdulist = pyfits.HDUList([output_hdu])
            output_hdr = output_hdulist[0].header
            output_hdr.update("DATE-OBS", str(rawfiledictionary["r" + 
                filenamestring0]["observation_time"]), 
                "observation date, UTC")
            output_hdr.update("HJULDATE", rawfiledictionary["r" + 
                filenamestring0]["observation_time_hjd"], 
                "observation date, Heliocentric Julian Date")
            output_hdr.update("RA", rawfiledictionary["r" + 
                filenamestring0]["observation_ra"], 
                "telescope pointing, low accuracy")
            output_hdr.update("DEC", rawfiledictionary["r" + 
                filenamestring0]["observation_dec"], 
                "telescope pointing, low accuracy")
            output_hdr.update("EPOCH", 2000.0, 
                "Epoch")
            output_hdr.update("EQUINOX", 2000.0, 
                "Equinox")
            output_hdr.update("SECPIX", 2.0, 
                "Plate scale, arcseconds per pixel")
            output_hdr.update("CROTA2", 180.0, 
                "Image Twist +AXIS2 W of N, J2000 (deg)")
            output_hdr.update("OBJECT", rawfiledictionary["r" + 
                filenamestring0]["pairitel_db_observation_id"], 
                "PAIRITEL observation ID")
            output_hdr.update("TRGTNAME", rawfiledictionary["r" + 
                filenamestring0]["target_name"], 
                "target name")
            output_hdr.update("FILTER", image_type[0])
            if image_type[2:] == "long":
                output_hdr.update("EXPTIME", 23.400, 
                    "exposure time, seconds")
            if image_type[2:] == "short":
                output_hdr.update("EXPTIME", 0.153, 
                    "exposure time, seconds")
            output_hdr.update("GAIN", 10.5, "Nominal gain (e-/ADU)")
            output_hdr.update("AIRMASS", rawfiledictionary["r" + 
                filenamestring0]["airmass"], "secant z")
            output_hdr.update("MEDSKARK", medskark_data, "Median skark value.")
            output_hdr.update("STRT_CPU", rawfiledictionary["r" + 
                filenamestring0]["strt_cpu"])
            output_hdr.update("STOP_CPU", rawfiledictionary["r" + 
                filenamestring2]["stop_cpu"])
            output_hdr.update("STRT0000", rawfiledictionary["r" + 
                filenamestring0]["strt_cpu"])
            output_hdr.update("STOP0000", rawfiledictionary["r" + 
                filenamestring0]["stop_cpu"])
            output_hdr.update("STRT0001", rawfiledictionary["r" + 
                filenamestring1]["strt_cpu"])
            output_hdr.update("STOP0001", rawfiledictionary["r" + 
                filenamestring1]["stop_cpu"])
            output_hdr.update("STRT0002", rawfiledictionary["r" + 
                filenamestring2]["strt_cpu"])
            output_hdr.update("STOP0002", rawfiledictionary["r" + 
                filenamestring2]["stop_cpu"])
            output_hdr.update("PI", rawfiledictionary["r" + 
                filenamestring0]["principal_investigator"], 
                "principal investigator")
            output_hdr.update("INSTRUME", "2MASS Survey cam")
            output_hdr.update("OBSERVAT", "Mt. Hopkins")
            output_hdr.update("TELESCOP", "1.3m PAIRITEL")
            output_hdr.update("REDUX-SW", "PAIRITEL Pipeline v3.5", 
                "Version of reduction software")
            output_hdr.update("JDETAREA", j_min_area, 
                "Min detection area for J object masking")
            output_hdr.update("HDETAREA", h_min_area, 
                "Min detection area for H object masking")
            output_hdr.update("KDETAREA", k_min_area, 
                "Min detection area for K object masking")
            output_hdr.update("JTHRESH", j_thresh, 
                "Min detection threshold for J object masking")
            output_hdr.update("HTHRESH", h_thresh, 
                "Min detection threshold for H object masking")
            output_hdr.update("KTHRESH", k_thresh, 
                "Min detection threshold for K object masking")            
                # Finally, just write out the image file.
            output_hdulist.writeto("./triplestacks/" + image_type + 
                "_triplestack" + filenamestring0[:-7] + ".fits")
system("mkdir triplestacks/cloudy")
system("mkdir triplestackweights/cloudy")
if do_cloud:
    for cloudy_triplestack in cloudy_triplestacks_list:
        if cloudy_triplestack[0].split("-p")[1] != "0.fits":
            image_string = (cloudy_triplestack[0].split("j_long_")[0] + "*" + 
                cloudy_triplestack[0].split("j_long_")[1])
            system("mv " + image_string + " triplestacks/cloudy")
            image_string = image_string.replace("triplestacks", 
                "triplestackweights")
            image_string = image_string.replace("*triplestack", 
                "*triplestackweightmap")
            system("mv " + image_string + " triplestackweights/cloudy")
# Run scamp on the source_imagefile. If it fails, can try anet.py.
#print "Now fitting WCS to triplestacks via scamp"
#source_weightfile = source_imagefile.replace("_triplestack", 
#    "_triplestackweightmap").replace("triplestacks", "triplestackweights")
#source_imagefile_hdulist = pyfits.open(source_imagefile)
#source_imagefile_header = source_imagefile_hdulist[0].header
#source_imagefile_CRVAL1 = float(source_imagefile_header["RA"])
#source_imagefile_CRVAL2 = float(source_imagefile_header["DEC"])
#source_imagefile_hdulist.close()
#system(sethead_bin + " " + source_imagefile + 
#    " WAT1_001=" + "wtype=tan" +
#    " WAT2_001=" + "wtype=tan" +
#    " RADECSYS=" + "FK5" +
#    " CTYPE1=" + "RA---TAN-SIP" +
#    " CTYPE2=" + "DEC---TAN-SIP" +
#    " CRPIX1=" + "128.5" +
#    " CRPIX2=" + "128.5" + 
#    " SECPIX1=" + "2.0" + 
#    " SECPIX2=" + "2.0" +
#    " CRVAL1=" + str(source_imagefile_CRVAL1) +
#    " CRVAL2=" + str(source_imagefile_CRVAL2) +
#    " CROTA1=" + "180.0" +
#    " CROTA2=" + "0.0" +
#    " CD1_1=" + "0.000550913348843" +
#    " CD1_2=" + "-7.56800716897e-06" +
#    " CD2_1=" + "-7.82421665954e-06" +
#    " CD2_2=" + "-0.000550553433267" +
#    " GAIN=10.5")
#system(sethead_bin + " " + source_weightfile + 
#    " WAT1_001=" + "wtype=tan" +
#    " WAT2_001=" + "wtype=tan" +
#    " RADECSYS=" + "FK5" +
#    " CTYPE1=" + "RA---TAN-SIP" +
#    " CTYPE2=" + "DEC---TAN-SIP" +
#    " CRPIX1=" + "128.5" +
#    " CRPIX2=" + "128.5" + 
#    " SECPIX1=" + "2.0" + 
#    " SECPIX2=" + "2.0" +
#    " CRVAL1=" + str(source_imagefile_CRVAL1) +
#    " CRVAL2=" + str(source_imagefile_CRVAL2) +
#    " CROTA1=" + "180.0" +
#    " CROTA2=" + "0.0" +
#    " CD1_1=" + "0.000550913348843" +
#    " CD1_2=" + "-7.56800716897e-06" +
#    " CD2_1=" + "-7.82421665954e-06" +
#    " CD2_2=" + "-0.000550553433267" +
#    " GAIN=10.5")
#print source_imagefile
#scamp_worked = scamp_call(source_imagefile, source_weightfile)
#print scamp_worked
# Run anet.py on the source_imagefile. Then copy the WCS header info
# into all the other triplestacks and reduced images.
print "Now fitting WCS to triplestacks via astrometry.net"
system(python_bin + " anet_triplestack.py " + source_imagefile)
source_imagefile_hdulist = pyfits.open(source_imagefile)
source_imagefile_header = source_imagefile_hdulist[0].header
try:
    source_imagefile_CRPIX1 = float(source_imagefile_header["CRPIX1"])
    astromery_dot_net_worked = True
except:
    source_imagefile_CRPIX1 = 128.5
    astromery_dot_net_worked = False
try:
    source_imagefile_CRPIX2 = float(source_imagefile_header["CRPIX2"])
except:
    source_imagefile_CRPIX2 = 128.5
try:
    source_imagefile_SECPIX1 = float(source_imagefile_header["SECPIX1"])
except:
    source_imagefile_SECPIX1 = 2.0
try:
    source_imagefile_SECPIX2 = float(source_imagefile_header["SECPIX2"])
except:
    source_imagefile_SECPIX2 = 2.0
source_imagefile_EPOCH = float(source_imagefile_header["EPOCH"])                   
source_imagefile_EQUINOX = float(source_imagefile_header["EQUINOX"])
try:
    source_imagefile_CRVAL1 = float(source_imagefile_header["CRVAL1"])
except:
    source_imagefile_CRVAL1 = float(source_imagefile_header["RA"])
try:
    source_imagefile_CRVAL2 = float(source_imagefile_header["CRVAL2"])
except:
    source_imagefile_CRVAL2 = float(source_imagefile_header["DEC"])
try:
    source_imagefile_CROTA1 = float(source_imagefile_header["CROTA1"])
except:
    source_imagefile_CROTA1 = 180.0
source_imagefile_CROTA2 = 0.0
try:
    source_imagefile_CD1_1 = float(source_imagefile_header["CD1_1"])
except:
    source_imagefile_CD1_1 = 0.00055
try:
    source_imagefile_CD1_2 = float(source_imagefile_header["CD1_2"])
except:
    source_imagefile_CD1_2 = 0.0
try:
    source_imagefile_CD2_1 = float(source_imagefile_header["CD2_1"])
except:
    source_imagefile_CD2_1 = 0.0
try:
    source_imagefile_CD2_2 = float(source_imagefile_header["CD2_2"])
except:
    source_imagefile_CD2_2 = -0.00055
source_imagefile_RA = source_imagefile_header["RA"]
source_imagefile_DEC = source_imagefile_header["DEC"]
try:
    source_imagefile_WAT1_001 = source_imagefile_header["WAT1_001"]
except:
    source_imagefile_WAT1_001 = "wtype=tan"
try:
    source_imagefile_WAT2_001 = source_imagefile_header["WAT2_001"]
except:
    source_imagefile_WAT2_001 = "wtype=tan"
try:
    source_imagefile_RADECSYS = source_imagefile_header["RADECSYS"]
except:
    source_imagefile_RADECSYS = "FK5"
try:
    source_imagefile_CTYPE1 = source_imagefile_header["CTYPE1"]
except:
    source_imagefile_CTYPE1 = "RA---TAN"
try:
    source_imagefile_CTYPE2 = source_imagefile_header["CTYPE2"]
except:
    source_imagefile_CTYPE2 = "DEC--TAN"
source_imagefile_hdulist.close()
print "Finished triplestacks, %f seconds required." % (time() - t1)
system(sethead_bin + " ./triplestacks/*.fits" + 
    " RA=" + str(source_imagefile_RA) +
    " DEC=" + str(source_imagefile_DEC) +
    " WAT1_001=" + source_imagefile_WAT1_001 +
    " WAT2_001=" + source_imagefile_WAT2_001 +
    " RADECSYS=" + source_imagefile_RADECSYS +
    " CTYPE1=" + source_imagefile_CTYPE1 +
    " CTYPE2=" + source_imagefile_CTYPE2 +
    " CRPIX1=" + str(source_imagefile_CRPIX1) +
    " CRPIX2=" + str(source_imagefile_CRPIX2) + 
    " SECPIX1=" + str(source_imagefile_SECPIX1) + 
    " SECPIX2=" + str(source_imagefile_SECPIX2) +
    " EPOCH=" + str(source_imagefile_EPOCH) +
    " EQUINOX=" + str(source_imagefile_EQUINOX) +
    " CRVAL1=" + str(source_imagefile_CRVAL1) +
    " CRVAL2=" + str(source_imagefile_CRVAL2) +
    " CROTA1=" + str(source_imagefile_CROTA1) +
    " CROTA2=" + str(source_imagefile_CROTA2) +
    " CD1_1=" + str(source_imagefile_CD1_1) +
    " CD1_2=" + str(source_imagefile_CD1_2) +
    " CD2_1=" + str(source_imagefile_CD2_1) +
    " CD2_2=" + str(source_imagefile_CD2_2) +
    " GAIN=10.5")
system(sethead_bin + " ./reduced/*.fits" + 
    " RA=" + str(source_imagefile_RA) +
    " DEC=" + str(source_imagefile_DEC) +
    " WAT1_001=" + source_imagefile_WAT1_001 +
    " WAT2_001=" + source_imagefile_WAT2_001 +
    " RADECSYS=" + source_imagefile_RADECSYS +
    " CTYPE1=" + source_imagefile_CTYPE1 +
    " CTYPE2=" + source_imagefile_CTYPE2 +
    " CRPIX1=" + str(source_imagefile_CRPIX1) +
    " CRPIX2=" + str(source_imagefile_CRPIX2) + 
    " SECPIX1=" + str(source_imagefile_SECPIX1) + 
    " SECPIX2=" + str(source_imagefile_SECPIX2) +
    " EPOCH=" + str(source_imagefile_EPOCH) +
    " EQUINOX=" + str(source_imagefile_EQUINOX) +
    " CRVAL1=" + str(source_imagefile_CRVAL1) +
    " CRVAL2=" + str(source_imagefile_CRVAL2) +
    " CROTA1=" + str(source_imagefile_CROTA1) +
    " CROTA2=" + str(source_imagefile_CROTA2) +
    " CD1_1=" + str(source_imagefile_CD1_1) +
    " CD1_2=" + str(source_imagefile_CD1_2) +
    " CD2_1=" + str(source_imagefile_CD2_1) +
    " CD2_2=" + str(source_imagefile_CD2_2) +
    " GAIN=10.5")
# Now, apply the calculated pixel offsets to the triplestacks and reduced 
# frames. The h and k bands have offset offsets because the detectors are not
# perfectingly aligned in the center of the lightpath of the camera. These 
# additional offsets account for this and ensure good WCS fits for the 
# triplestacks and reduced images in all bands.
alignment_file = file("final_alignments_j.txt", "r")
sethead_file = file("sethead_commands.txt", "w")
for line in alignment_file:
    xoffset = float(line.split()[0])
    yoffset = float(line.split()[1])
    image_file = line.split()[2]
    sethead_file.write(sethead_bin + " " + image_file.replace("j_long", "j_*") + 
        " CRPIX1=" + str(source_imagefile_CRPIX1 - xoffset) + " CRPIX2=" + 
        str(source_imagefile_CRPIX2 + yoffset) + 
        " \n")
    sethead_file.write(sethead_bin + " " + 
        image_file.replace("j_long", "j_*").replace("triplestacks", 
        "reduced").replace("triplestack", "reduced").replace(".fits", "-*.fits")
        + " CRPIX1=" + str(source_imagefile_CRPIX1 - xoffset) + " CRPIX2=" + 
        str(source_imagefile_CRPIX2 + yoffset) + 
        " \n")
    sethead_file.write(sethead_bin + " " + image_file.replace("j_long", "h_*") + 
        " CRPIX1=" + str(source_imagefile_CRPIX1 - xoffset + 1.246074) + 
        " CRPIX2=" + str(source_imagefile_CRPIX2 + yoffset + 0.0484713) + 
        " \n")
    sethead_file.write(sethead_bin + " " + 
        image_file.replace("j_long", "h_*").replace("triplestacks", 
        "reduced").replace("triplestack", "reduced").replace(".fits", "-*.fits")
        + " CRPIX1=" + str(source_imagefile_CRPIX1 - xoffset + 1.246074) + 
        " CRPIX2=" + str(source_imagefile_CRPIX2 + yoffset + 0.0484713) + 
        " \n")
    sethead_file.write(sethead_bin + " " + image_file.replace("j_long", "k_*") + 
        " CRPIX1=" + str(source_imagefile_CRPIX1 - xoffset + 0.2991059) + 
        " CRPIX2=" + str(source_imagefile_CRPIX2 + yoffset - 3.039885) + 
        " \n")
    sethead_file.write(sethead_bin + " " + 
        image_file.replace("j_long", "k_*").replace("triplestacks", 
        "reduced").replace("triplestack", "reduced").replace(".fits", "-*.fits")
        + " CRPIX1=" + str(source_imagefile_CRPIX1 - xoffset + 0.2991059) + 
        " CRPIX2=" + str(source_imagefile_CRPIX2 + yoffset - 3.039885) + 
        " \n")
alignment_file.close()
sethead_file.close()
system("chmod +x sethead_commands.txt")
system("./sethead_commands.txt")
system("mkdir triplestacks/bad_alignments")
system("mkdir triplestackweights/bad_alignments")
system("mkdir reduced/bad_alignments")
num_bad_alignments = 0
alignmets_evaluation_file = file("alignments_evaluation.txt", "w")
for n in range(len(final_alignments_list_sorted)):
    if final_alignments_list_sorted[n][2] == pointing_alignments_list[n][2]:
        pixel_diff = sqrt((float(final_alignments_list_sorted[n][0]) - 
            float(pointing_alignments_list[n][0]))**2 + 
            (float(final_alignments_list_sorted[n][1]) - 
            float(pointing_alignments_list[n][1]))**2)
        alignmets_evaluation_file.write(str(pixel_diff) + " " + 
            final_alignments_list_sorted[n][2] + "\n")
        if pixel_diff > 15 or pixel_diff < 0:
            image_string = str(final_alignments_list_sorted[n][2])
            image_string = image_string.replace("j_long", "*")
            if (str(final_alignments_list_sorted[n][2]).split("-p")[1] !=
                "0.fits"):
                system("mv " + image_string + " triplestacks/bad_alignments")
            system("mv " + 
                image_string.replace("triplestacks", 
                "reduced").replace("triplestack", "reduced")[:-5] + 
                "* reduced/bad_alignments")
            if (str(final_alignments_list_sorted[n][2]).split("-p")[1] !=
                "0.fits"):
                image_string = image_string.replace("triplestacks/", 
                    "triplestackweights/")
                image_string = image_string.replace("_triplestack", 
                    "_triplestackweightmap")
                system("mv " + image_string + 
                    " triplestackweights/bad_alignments")
            num_bad_alignments = num_bad_alignments + 1
alignmets_evaluation_file.close()
if num_bad_alignments > 0:
    print str(num_bad_alignments) + " bad alignments detected, rejecting."
if do_short:
    j_short_triplestacks = file("j_short_triplestacks.txt", "w")
    h_short_triplestacks = file("h_short_triplestacks.txt", "w")
    k_short_triplestacks = file("k_short_triplestacks.txt", "w")
j_long_triplestacks = file("j_long_triplestacks.txt", "w")
h_long_triplestacks = file("h_long_triplestacks.txt", "w")
k_long_triplestacks = file("k_long_triplestacks.txt", "w")
j_long_triplestacks_list = []
triplestackfilelistfull = listdir("./triplestacks")
for filename in triplestackfilelistfull:
    if filename[-5:] == ".fits" and filename.split("-")[-1] != "p0.fits":
        if do_short:
            if filename[:7] == "j_short":
                j_short_triplestacks.write("triplestacks/" + filename + "\n")
            if filename[:7] == "h_short":
                h_short_triplestacks.write("triplestacks/" + filename + "\n")
            if filename[:7] == "k_short":
                k_short_triplestacks.write("triplestacks/" + filename + "\n")
        if filename[:6] == "j_long":
            j_long_triplestacks.write("triplestacks/" + filename + "\n")
            j_long_triplestacks_list.append(filename)
        if filename[:6] == "h_long":
            h_long_triplestacks.write("triplestacks/" + filename + "\n")
        if filename[:6] == "k_long":
            k_long_triplestacks.write("triplestacks/" + filename + "\n")
if do_short:
    j_short_triplestacks.close()
    h_short_triplestacks.close()
    k_short_triplestacks.close()
j_long_triplestacks.close()
h_long_triplestacks.close()
k_long_triplestacks.close()
if do_short:
    j_short_triplestackweights = file("j_short_triplestackweights.txt", "w")
    h_short_triplestackweights = file("h_short_triplestackweights.txt", "w")
    k_short_triplestackweights = file("k_short_triplestackweights.txt", "w")
j_long_triplestackweights = file("j_long_triplestackweights.txt", "w")
h_long_triplestackweights = file("h_long_triplestackweights.txt", "w")
k_long_triplestackweights = file("k_long_triplestackweights.txt", "w")
triplestackfilelistfull = listdir("./triplestackweights")
for filename in triplestackfilelistfull:
    if filename[-5:] == ".fits" and filename.split("-")[-1] != "p0.fits":
        if do_short:
            if filename[:7] == "j_short":
                j_short_triplestackweights.write("triplestackweights/" + 
                    filename + "\n")
            if filename[:7] == "h_short":
                h_short_triplestackweights.write("triplestackweights/" + 
                    filename + "\n")
            if filename[:7] == "k_short":
                k_short_triplestackweights.write("triplestackweights/" + 
                    filename + "\n")
        if filename[:6] == "j_long":
            j_long_triplestackweights.write("triplestackweights/" + 
                filename + "\n")
        if filename[:6] == "h_long":
            h_long_triplestackweights.write("triplestackweights/" + 
                filename + "\n")
        if filename[:6] == "k_long":
            k_long_triplestackweights.write("triplestackweights/" + 
                filename + "\n")
if do_short:
    j_short_triplestackweights.close()
    h_short_triplestackweights.close()
    k_short_triplestackweights.close()
j_long_triplestackweights.close()
h_long_triplestackweights.close()
k_long_triplestackweights.close()
print "Finished aligning images, %f seconds required." % (time() - t1)
# Run the mosaicing.
print "Now running image coaddition with SWarp."
t1 = time()
# So far, we've reduce the p0 triplestacks as normal (although, we didn't 
# calculate or apply offsets). We don't want to use the p0 triplestacks to make
# the coadd mosaics becuase the p0-0 frames are usually dominated by read noise.
# This is likely due to the long telescope slew from the last observation. So,
# move the p0 triplestacks to a subdirectory before running SWarp.
system("mkdir ./triplestacks/p0_triplestacks")
system("mv ./triplestacks/*-p0.fits ./triplestacks/p0_triplestacks")
system("mkdir ./triplestackweights/p0_triplestackweights")
system("mv ./triplestackweights/*-p0.fits " + 
    "./triplestackweights/p0_triplestackweights")       
# Make list of swarp_commands.
if do_short:
    if obs_string.split(".")[0] == "SN" or obs_string.split(".")[0] == "SNTEMP":
        swarp_commands = [
            swarp_bin + " @j_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @j_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME j_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME j_long_" + obs_string + "_coadd.weight.fits " + 
            "-IMAGE_SIZE 800,800",
            swarp_bin + " @j_short_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @j_short_triplestackweights.txt " + 
            "-IMAGEOUT_NAME j_short_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME j_short_" + obs_string + "_coadd.weight.fits " + 
            "-IMAGE_SIZE 800,800",
            swarp_bin + " @h_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @h_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME h_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME h_long_" + obs_string + "_coadd.weight.fits " + 
            "-IMAGE_SIZE 800,800",
            swarp_bin + " @h_short_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @h_short_triplestackweights.txt " + 
            "-IMAGEOUT_NAME h_short_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME h_short_" + obs_string + "_coadd.weight.fits " + 
            "-IMAGE_SIZE 800,800",
            swarp_bin + " @k_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @k_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME k_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME k_long_" + obs_string + "_coadd.weight.fits " + 
            "-IMAGE_SIZE 800,800",
            swarp_bin + " @k_short_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @k_short_triplestackweights.txt " + 
            "-IMAGEOUT_NAME k_short_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME k_short_" + obs_string + "_coadd.weight.fits " + 
            "-IMAGE_SIZE 800,800"]
    else:
        swarp_commands = [
            swarp_bin + " @j_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @j_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME j_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME j_long_" + obs_string + "_coadd.weight.fits",
            swarp_bin + " @j_short_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @j_short_triplestackweights.txt " + 
            "-IMAGEOUT_NAME j_short_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME j_short_" + obs_string + "_coadd.weight.fits",
            swarp_bin + " @h_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @h_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME h_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME h_long_" + obs_string + "_coadd.weight.fits",
            swarp_bin + " @h_short_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @h_short_triplestackweights.txt " + 
            "-IMAGEOUT_NAME h_short_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME h_short_" + obs_string + "_coadd.weight.fits",
            swarp_bin + " @k_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @k_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME k_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME k_long_" + obs_string + "_coadd.weight.fits",
            swarp_bin + " @k_short_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @k_short_triplestackweights.txt " + 
            "-IMAGEOUT_NAME k_short_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME k_short_" + obs_string + "_coadd.weight.fits"]

if not do_short:
    if obs_string.split(".")[0] == "SN" or obs_string.split(".")[0] == "SNTEMP":
        swarp_commands = [
            swarp_bin + " @j_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @j_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME j_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME j_long_" + obs_string + "_coadd.weight.fits " + 
            "-IMAGE_SIZE 800,800",
            swarp_bin + " @h_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @h_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME h_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME h_long_" + obs_string + "_coadd.weight.fits " + 
            "-IMAGE_SIZE 800,800",
            swarp_bin + " @k_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @k_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME k_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME k_long_" + obs_string + "_coadd.weight.fits " + 
            "-IMAGE_SIZE 800,800"]
    else:
        swarp_commands = [
            swarp_bin + " @j_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @j_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME j_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME j_long_" + obs_string + "_coadd.weight.fits",
            swarp_bin + " @h_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @h_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME h_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME h_long_" + obs_string + "_coadd.weight.fits",
            swarp_bin + " @k_long_triplestacks.txt " + 
            "-c pairitel_redux.swarp " + 
            "-WEIGHT_IMAGE @k_long_triplestackweights.txt " + 
            "-IMAGEOUT_NAME k_long_" + obs_string + "_coadd.fits " + 
            "-WEIGHTOUT_NAME k_long_" + obs_string + "_coadd.weight.fits"]
# Run the mosaicing with parallel processing.
p = Pool(numprocessors)
result = p.map_async(run_swarp, swarp_commands)
poolresult = result.get()
# We insert the STRT_CPU of the first triplestack and the STOP_CPU of 
# the last triplestack used to make the mosaic.
j_long_triplestacks_list_pnum = []
for j_long_triplestack in j_long_triplestacks_list:
    pnum = int(j_long_triplestack.split("-p")[1].split(".")[0])
    j_long_triplestacks_list_pnum.append([j_long_triplestack, pnum])
j_long_triplestacks_list_pnum_sorted = sorted(j_long_triplestacks_list_pnum,
    key=operator.itemgetter(1))
j_long_triplestacks_list = []
for entry in j_long_triplestacks_list_pnum_sorted:
    j_long_triplestacks_list.append(entry[0])
trip_hdulist = pyfits.open("./triplestacks/" + j_long_triplestacks_list[0])
trip_header = trip_hdulist[0].header
strt_cpu = str(trip_header["STRT_CPU"])
trip_hdulist.close()
trip_hdulist = pyfits.open("./triplestacks/" + j_long_triplestacks_list[-1])
trip_header = trip_hdulist[0].header
stop_cpu = str(trip_header["STOP_CPU"])
trip_hdulist.close()
system(sethead_bin + " *" + obs_string + "_coadd*.fits " + 
            "STRT_CPU='%s' STOP_CPU='%s'" % (strt_cpu, stop_cpu))
# And, finally, we insert the STRT_CPU and STOP_CPU for each rawfile which was
# used in the mosaic.
for triplestackfile in j_long_triplestacks_list:
    if ((triplestackfile[-5:] == ".fits") and 
        (triplestackfile[:6] == "j_long") and
        (triplestackfile.split("-p")[1].split(".")[0] != "0")):
        p_num = int(triplestackfile.split("-p")[1].split(".")[0])
        num_0 = ("%4.f" % ((p_num - 1) * 3)).replace(" ", "0")
        num_1 = ("%4.f" % ((p_num - 1) * 3 + 1)).replace(" ", "0")
        num_2 = ("%4.f" % ((p_num - 1) * 3 + 2)).replace(" ", "0")
        strt_0 = "STRT" + num_0
        strt_1 = "STRT" + num_1
        strt_2 = "STRT" + num_2
        stop_0 = "STOP" + num_0
        stop_1 = "STOP" + num_1
        stop_2 = "STOP" + num_2
        trip_hdulist = pyfits.open("./triplestacks/" + triplestackfile)
        trip_header = trip_hdulist[0].header
        strt_0_val = str(trip_header["STRT0000"])
        strt_1_val = str(trip_header["STRT0001"])
        strt_2_val = str(trip_header["STRT0002"])
        stop_0_val = str(trip_header["STOP0000"])
        stop_1_val = str(trip_header["STOP0001"])
        stop_2_val = str(trip_header["STOP0002"])
        trip_hdulist.close()
        system(sethead_bin + " *" + obs_string + "_coadd*.fits "
            "%s='%s' %s='%s' %s='%s' %s='%s' %s='%s' %s='%s'" % (strt_0, 
            strt_0_val, stop_0, stop_0_val, strt_1, strt_1_val, stop_1, 
            stop_1_val, strt_2, strt_2_val, stop_2, stop_2_val))

j_long_meds = []
h_long_meds = []
k_long_meds = []
if do_short:
    j_short_meds = []
    h_short_meds = []
    k_short_meds = []

for triplestackfile in j_long_triplestacks_list:
    if ((triplestackfile[-5:] == ".fits") and 
        (triplestackfile[:6] == "j_long") and
        (triplestackfile.split("-p")[1].split(".")[0] != "0")):
        trip_hdulist = pyfits.open("./triplestacks/" + triplestackfile)
        trip_header = trip_hdulist[0].header
        j_long_meds.append(float(trip_header["MEDSKARK"]))
        trip_hdulist.close()
        trip_hdulist = pyfits.open("./triplestacks/" + triplestackfile.replace("j_long", "h_long"))
        trip_header = trip_hdulist[0].header
        h_long_meds.append(float(trip_header["MEDSKARK"]))
        trip_hdulist.close()
        trip_hdulist = pyfits.open("./triplestacks/" + triplestackfile.replace("j_long", "k_long"))
        trip_header = trip_hdulist[0].header
        k_long_meds.append(float(trip_header["MEDSKARK"]))
        trip_hdulist.close()
        if do_short:
            trip_hdulist = pyfits.open("./triplestacks/" + triplestackfile.replace("j_long", "j_short"))
            trip_header = trip_hdulist[0].header
            j_short_meds.append(float(trip_header["MEDSKARK"]))
            trip_hdulist.close()
            trip_hdulist = pyfits.open("./triplestacks/" + triplestackfile.replace("j_long", "h_short"))
            trip_header = trip_hdulist[0].header
            h_short_meds.append(float(trip_header["MEDSKARK"]))
            trip_hdulist.close()
            trip_hdulist = pyfits.open("./triplestacks/" + triplestackfile.replace("j_long", "k_short"))
            trip_header = trip_hdulist[0].header
            k_short_meds.append(float(trip_header["MEDSKARK"]))
            trip_hdulist.close()
            
j_long_mosaic_medskark = median(j_long_meds)
system(sethead_bin + " j_long_" + obs_string + "_coadd.fits MEDSKARK=" + str(j_long_mosaic_medskark))
h_long_mosaic_medskark = median(h_long_meds)
system(sethead_bin + " h_long_" + obs_string + "_coadd.fits MEDSKARK=" + str(h_long_mosaic_medskark))
k_long_mosaic_medskark = median(k_long_meds)
system(sethead_bin + " k_long_" + obs_string + "_coadd.fits MEDSKARK=" + str(k_long_mosaic_medskark))
if do_short:
    j_short_mosaic_medskark = median(j_short_meds)
    system(sethead_bin + " j_short_" + obs_string + "_coadd.fits MEDSKARK=" + str(j_short_mosaic_medskark))
    h_short_mosaic_medskark = median(h_short_meds)
    system(sethead_bin + " h_short_" + obs_string + "_coadd.fits MEDSKARK=" + str(h_short_mosaic_medskark))
    k_short_mosaic_medskark = median(k_short_meds)
    system(sethead_bin + " k_short_" + obs_string + "_coadd.fits MEDSKARK=" + str(k_short_mosaic_medskark))

print "Finished image coaddition, %f seconds required." % (time() - t1)
# Run the wcs fitting.
print "Now wcs fitting coadd images."
t1 = time()
scamp_worked = False
scamp_worked = scamp_call("j_long_" + obs_string + "_coadd.fits", 
    "j_long_" + obs_string + "_coadd.weight.fits")
if scamp_worked:
    system(python_bin + " missfits2.py j_long_" + obs_string + "_coadd.head " + 
        "j_long_" + obs_string + "_coadd.weight.fits")
    system(python_bin + " missfits2.py h_long_" + obs_string + "_coadd.head " + 
        "h_long_" + obs_string + "_coadd.weight.fits")
    system(python_bin + " missfits2.py k_long_" + obs_string + "_coadd.head " + 
        "k_long_" + obs_string + "_coadd.weight.fits")
    if do_short:
        system(python_bin + " missfits2.py j_long_" + obs_string + "_coadd.head " + 
            "j_short_" + obs_string + "_coadd.fits")
        system(python_bin + " missfits2.py j_long_" + obs_string + "_coadd.head " + 
            "j_short_" + obs_string + "_coadd.weight.fits")
        system(python_bin + " missfits2.py h_long_" + obs_string + "_coadd.head " + 
            "h_short_" + obs_string + "_coadd.fits")
        system(python_bin + " missfits2.py h_long_" + obs_string + "_coadd.head " + 
            "h_short_" + obs_string + "_coadd.weight.fits")
        system(python_bin + " missfits2.py k_long_" + obs_string + "_coadd.head " + 
            "k_short_" + obs_string + "_coadd.fits")
        system(python_bin + " missfits2.py k_long_" + obs_string + "_coadd.head " + 
            "k_short_" + obs_string + "_coadd.weight.fits")
system("rm *.head")
system("rm *.ldac")
if not scamp_worked:
    system(python_bin + " anet.py *_long_" + obs_string + "_coadd.fits")
    j_long_hdulist = pyfits.open("j_long_" + obs_string + "_coadd.fits", 
        "readonly")
    j_long_weights_hdulist = pyfits.open("j_long_" + obs_string + 
        "_coadd.weight.fits", "update")
    if do_short:
        j_short_hdulist = pyfits.open("j_short_" + obs_string + "_coadd.fits", 
            "update")
        j_short_weights_hdulist = pyfits.open("j_short_" + obs_string + 
            "_coadd.weight.fits", "update")
    tmp0 = j_long_hdulist[0].header
    try:
        j_long_an_jobid = tmp0["AN_JOBID"]
    except:
        j_long_an_jobid = False
    del tmp0["SIMPLE"]
    del tmp0["BITPIX"]
    del tmp0["NAXIS"]
    j_long_hdulist.close()
    h_long_hdulist = pyfits.open("h_long_" + obs_string + "_coadd.fits", 
        "readonly")
    h_long_weights_hdulist = pyfits.open("h_long_" + obs_string + 
        "_coadd.weight.fits", "update")
    if do_short:
        h_short_hdulist = pyfits.open("h_short_" + obs_string + "_coadd.fits", 
            "update")
        h_short_weights_hdulist = pyfits.open("h_short_" + obs_string + 
            "_coadd.weight.fits", "update")
    tmp1 = h_long_hdulist[0].header
    try:
        h_long_an_jobid = tmp1["AN_JOBID"]
    except:
        h_long_an_jobid = False
    del tmp1["SIMPLE"]
    del tmp1["BITPIX"]
    del tmp1["NAXIS"]
    h_long_hdulist.close()
    k_long_hdulist = pyfits.open("k_long_" + obs_string + "_coadd.fits", 
        "readonly")
    k_long_weights_hdulist = pyfits.open("k_long_" + obs_string + 
        "_coadd.weight.fits", "update")
    if do_short:
        k_short_hdulist = pyfits.open("k_short_" + obs_string + "_coadd.fits", 
            "update")
        k_short_weights_hdulist = pyfits.open("k_short_" + obs_string + 
            "_coadd.weight.fits", "update")
    tmp2 = k_long_hdulist[0].header
    try:
        k_long_an_jobid = tmp2["AN_JOBID"]
    except:
        k_long_an_jobid = False
    del tmp2["SIMPLE"]
    del tmp2["BITPIX"]
    del tmp2["NAXIS"]
    k_long_hdulist.close()
    if j_long_an_jobid:
        for c in tmp0.ascardlist():
            j_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
            if do_short:
                j_short_hdulist[0].header.update(c.key,c.value,c.comment)
                j_short_weights_hdulist[0].header.update(c.key,c.value,c.comment)
    if h_long_an_jobid:
        for c in tmp1.ascardlist():
            h_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
            if do_short:
                h_short_hdulist[0].header.update(c.key,c.value,c.comment)
                h_short_weights_hdulist[0].header.update(c.key,c.value,c.comment)
    if k_long_an_jobid:
        for c in tmp2.ascardlist():
            k_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
            if do_short:
                k_short_hdulist[0].header.update(c.key,c.value,c.comment)
                k_short_weights_hdulist[0].header.update(c.key,c.value,c.comment)
    if j_long_an_jobid and not h_long_an_jobid:
        h_long_hdulist = pyfits.open("h_long_" + obs_string + "_coadd.fits", 
        "update")
        for c in tmp0.ascardlist():
            if c.key == "FILTER":
                continue
            if c.key == "CRPIX1":
                h_long_hdulist[0].header.update(c.key,float(c.value) - 
                    2.492148,c.comment)
                h_long_weights_hdulist[0].header.update(c.key,float(c.value) - 
                    2.492148,c.comment)
                if do_short:
                    h_short_hdulist[0].header.update(c.key,float(c.value) - 
                        2.492148,c.comment)
                    h_short_weights_hdulist[0].header.update(c.key,float(c.value) - 
                        2.492148,c.comment)
            elif c.key == "CRPIX2":
                h_long_hdulist[0].header.update(c.key,float(c.value) - 
                    0.0969426,c.comment)
                h_long_weights_hdulist[0].header.update(c.key,float(c.value) - 
                    0.0969426,c.comment)
                if do_short:
                    h_short_hdulist[0].header.update(c.key,float(c.value) - 
                        0.0969426,c.comment)
                    h_short_weights_hdulist[0].header.update(c.key,float(c.value) - 
                        0.0969426,c.comment)
            else:
                h_long_hdulist[0].header.update(c.key,c.value,c.comment)
                h_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
                if do_short:
                    h_short_hdulist[0].header.update(c.key,c.value,c.comment)
                    h_short_weights_hdulist[0].header.update(c.key,c.value,
                        c.comment)
        h_long_hdulist.verify("silentfix")
        h_long_hdulist.close(output_verify='warn')  
    if j_long_an_jobid and not k_long_an_jobid:
        k_long_hdulist = pyfits.open("k_long_" + obs_string + "_coadd.fits", 
        "update")
        for c in tmp0.ascardlist():
            if c.key == "FILTER":
                continue
            elif c.key == "CRPIX1":
                k_long_hdulist[0].header.update(c.key,float(c.value) - 
                    0.5982118,c.comment)
                k_long_weights_hdulist[0].header.update(c.key,float(c.value) - 
                    0.5982118,c.comment)
                if do_short:
                    k_short_hdulist[0].header.update(c.key,float(c.value) - 
                        0.5982118,c.comment)
                    k_short_weights_hdulist[0].header.update(c.key,float(c.value) - 
                        0.5982118,c.comment)
            elif c.key == "CRPIX2":
                k_long_hdulist[0].header.update(c.key,float(c.value) + 
                    6.07977,c.comment)
                k_long_weights_hdulist[0].header.update(c.key,float(c.value) + 
                    6.07977,c.comment)
                if do_short:
                    k_short_hdulist[0].header.update(c.key,float(c.value) + 
                        6.07977,c.comment)
                    k_short_weights_hdulist[0].header.update(c.key,float(c.value) + 
                        6.07977,c.comment)
            elif c.key == "CD1_2":
                k_long_hdulist[0].header.update(c.key, 3.87797365654E-06, 
                    c.comment)
                k_long_weights_hdulist[0].header.update(c.key, 3.87797365654E-06, 
                    c.comment)
                if do_short:
                    k_short_hdulist[0].header.update(c.key, 3.87797365654E-06, 
                        c.comment)
                    k_short_weights_hdulist[0].header.update(c.key, 
                        3.87797365654E-06, c.comment)
            elif c.key == "CD2_1":
                k_long_hdulist[0].header.update(c.key, 4.22330716354E-06, 
                    c.comment)
                k_long_weights_hdulist[0].header.update(c.key, 4.22330716354E-06, 
                    c.comment)
                if do_short:
                    k_short_hdulist[0].header.update(c.key, 4.22330716354E-06, 
                        c.comment)
                    k_short_weights_hdulist[0].header.update(c.key, 
                        4.22330716354E-06, c.comment)
            else:
                k_long_hdulist[0].header.update(c.key,c.value,c.comment)
                k_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
                if do_short:
                    k_short_hdulist[0].header.update(c.key,c.value,c.comment)
                    k_short_weights_hdulist[0].header.update(c.key,c.value,
                        c.comment)
        k_long_hdulist.verify("silentfix")
        k_long_hdulist.close(output_verify='warn')
    j_long_weights_hdulist.verify("silentfix")
    j_long_weights_hdulist.close(output_verify='warn')
    if do_short:
        j_short_hdulist.verify("silentfix")
        j_short_hdulist.close(output_verify='warn')
        j_short_weights_hdulist.verify("silentfix")
        j_short_weights_hdulist.close(output_verify='warn')
    h_long_weights_hdulist.verify("silentfix")
    h_long_weights_hdulist.close(output_verify='warn')
    if do_short:
        h_short_hdulist.verify("silentfix")
        h_short_hdulist.close(output_verify='warn')
        h_short_weights_hdulist.verify("silentfix")
        h_short_weights_hdulist.close(output_verify='warn')
    k_long_weights_hdulist.verify("silentfix")
    k_long_weights_hdulist.close(output_verify='warn')
    if do_short:
        k_short_hdulist.verify("silentfix")
        k_short_hdulist.close(output_verify='warn')
        k_short_weights_hdulist.verify("silentfix")
        k_short_weights_hdulist.close(output_verify='warn')
# Make the triplestackfilelist.
triplestackfilelist = []
j_long_triplestacks_file = file("j_long_triplestacks.txt", "r")
for line in j_long_triplestacks_file:
    triplestackfilelist.append(line)
# We want to propagate the correct exposure time.
num_dither_positions = len(triplestackfilelist)
mosaic_long_exptime = num_dither_positions * 23.400
mosaic_short_exptime = num_dither_positions * 0.153
# Set the exposure time in the mosaics.
system(sethead_bin + " j_long_" + obs_string + "_coadd.fits EXPTIME=" + 
    str(mosaic_long_exptime))
if do_short:
    system(sethead_bin + " j_short_" + obs_string + "_coadd.fits EXPTIME=" + 
        str(mosaic_short_exptime))
system(sethead_bin + " h_long_" + obs_string + "_coadd.fits EXPTIME=" + 
    str(mosaic_long_exptime))
if do_short:
    system(sethead_bin + " h_short_" + obs_string + "_coadd.fits EXPTIME=" + 
        str(mosaic_short_exptime))
system(sethead_bin + " k_long_" + obs_string + "_coadd.fits EXPTIME=" + 
    str(mosaic_long_exptime))
if do_short:
    system(sethead_bin + " k_short_" + obs_string + "_coadd.fits EXPTIME=" + 
        str(mosaic_short_exptime))
        
# Set the RA and DEC header keywords in the final mosaics.
mosaic_hdulist = pyfits.open("j_long_" + obs_string + "_coadd.fits")
mosaic_header = mosaic_hdulist[0].header
crval1 = float(mosaic_header["CRVAL1"])
crval2 = float(mosaic_header["CRVAL2"])
crpix1 = float(mosaic_header["CRPIX1"])
crpix2 = float(mosaic_header["CRPIX2"])
mosaic_hdulist.close()
mosaic_dec = crval2 + (256.0 - crpix2)/3600
mosaic_ra = crval1 + (crpix1 - 256.0)/3600/cos(mosaic_dec*0.0174532925)
mosaic_dec = str(ephem.degrees(mosaic_dec * 0.0174532925))
mosaic_ra = str(ephem.hours(mosaic_ra * 0.0174532925))
system(sethead_bin + " j_long_" + obs_string + "_coadd.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " j_long_" + obs_string + "_coadd.weight.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " j_long_" + obs_string + "_coadd.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
system(sethead_bin + " j_long_" + obs_string + "_coadd.weight.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
if do_short:
    system(sethead_bin + " j_short_" + obs_string + "_coadd.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " j_short_" + obs_string + "_coadd.weight.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " j_short_" + obs_string + "_coadd.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
    system(sethead_bin + " j_short_" + obs_string + "_coadd.weight.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')

mosaic_hdulist = pyfits.open("h_long_" + obs_string + "_coadd.fits")
mosaic_header = mosaic_hdulist[0].header
crval1 = float(mosaic_header["CRVAL1"])
crval2 = float(mosaic_header["CRVAL2"])
crpix1 = float(mosaic_header["CRPIX1"])
crpix2 = float(mosaic_header["CRPIX2"])
mosaic_hdulist.close()
mosaic_dec = crval2 + (256.0 - crpix2)/3600
mosaic_ra = crval1 + (crpix1 - 256.0)/3600/cos(mosaic_dec*0.0174532925)
mosaic_dec = str(ephem.degrees(mosaic_dec * 0.0174532925))
mosaic_ra = str(ephem.hours(mosaic_ra * 0.0174532925))
system(sethead_bin + " h_long_" + obs_string + "_coadd.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " h_long_" + obs_string + "_coadd.weight.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " h_long_" + obs_string + "_coadd.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
system(sethead_bin + " h_long_" + obs_string + "_coadd.weight.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
if do_short:
    system(sethead_bin + " h_short_" + obs_string + "_coadd.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " h_short_" + obs_string + "_coadd.weight.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " h_short_" + obs_string + "_coadd.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
    system(sethead_bin + " h_short_" + obs_string + "_coadd.weight.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
        
mosaic_hdulist = pyfits.open("k_long_" + obs_string + "_coadd.fits")
mosaic_header = mosaic_hdulist[0].header
crval1 = float(mosaic_header["CRVAL1"])
crval2 = float(mosaic_header["CRVAL2"])
crpix1 = float(mosaic_header["CRPIX1"])
crpix2 = float(mosaic_header["CRPIX2"])
mosaic_hdulist.close()
mosaic_dec = crval2 + (256.0 - crpix2)/3600
mosaic_ra = crval1 + (crpix1 - 256.0)/3600/cos(mosaic_dec*0.0174532925)
mosaic_dec = str(ephem.degrees(mosaic_dec * 0.0174532925))
mosaic_ra = str(ephem.hours(mosaic_ra * 0.0174532925))
system(sethead_bin + " k_long_" + obs_string + "_coadd.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " k_long_" + obs_string + "_coadd.weight.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " k_long_" + obs_string + "_coadd.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
system(sethead_bin + " k_long_" + obs_string + "_coadd.weight.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
if do_short:
    system(sethead_bin + " k_short_" + obs_string + "_coadd.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " k_short_" + obs_string + "_coadd.weight.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " k_short_" + obs_string + "_coadd.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
    system(sethead_bin + " k_short_" + obs_string + "_coadd.weight.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
        
# Set the WCSFIT header keyword.
mosaic_hdulist = pyfits.open("j_long_" + obs_string + "_coadd.fits")
mosaic_header = mosaic_hdulist[0].header
try:
    scamp_solved = float(mosaic_header["PV1_1"])
    if scamp_solved > 0:
        scamp_solved = True
except:
    scamp_solved = False
try:
    an_solved = str(mosaic_header["A_ORDER"])
    if an_solved > 0:
        an_solved = True
except:
    an_solved = False
mosaic_hdulist.close()
if scamp_solved or an_solved:
    system(sethead_bin + " j_long_" + obs_string + "_coadd.fits WCSFIT=True")
    system(sethead_bin + " j_long_" + obs_string + "_coadd.weight.fits WCSFIT=True")
    if do_short:
        system(sethead_bin + " j_short_" + obs_string + "_coadd.fits WCSFIT=True")
        system(sethead_bin + " j_short_" + obs_string + "_coadd.weight.fits WCSFIT=True")
else:
    system(sethead_bin + " j_long_" + obs_string + "_coadd.fits WCSFIT=False")
    try:
        system(python_bin + " autoastrometry.py j_long_" + obs_string + "_coadd.fits -px 1.00 -inv -c ub2")
        j_long_hdulist = pyfits.open("j_long_" + obs_string + "_coadd.fits", "readonly")
        j_long_weights_hdulist = pyfits.open("j_long_" + obs_string + 
            "_coadd.weight.fits", "update")
        if do_short:
            j_short_hdulist = pyfits.open("j_short_" + obs_string + "_coadd.fits", 
                "update")
            j_short_weights_hdulist = pyfits.open("j_short_" + obs_string + 
                "_coadd.weight.fits", "update")
        tmp0 = j_long_hdulist[0].header
        del tmp0["SIMPLE"]
        del tmp0["BITPIX"]
        del tmp0["NAXIS"]
        j_long_hdulist.close()
        for c in tmp0.ascardlist():
            j_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
            if do_short:
                j_short_hdulist[0].header.update(c.key,c.value,c.comment)
                j_short_weights_hdulist[0].header.update(c.key,c.value,c.comment)
        j_long_weights_hdulist.verify("silentfix")
        j_long_weights_hdulist.close(output_verify='warn')
        if do_short:
            j_short_hdulist.verify("silentfix")
            j_short_hdulist.close(output_verify='warn')
            j_short_weights_hdulist.verify("silentfix")
            j_short_weights_hdulist.close(output_verify='warn')
    except:
        print "Failed J WCS fit."
        
mosaic_hdulist = pyfits.open("h_long_" + obs_string + "_coadd.fits")
mosaic_header = mosaic_hdulist[0].header
try:
    scamp_solved = float(mosaic_header["PV1_1"])
    if scamp_solved > 0:
        scamp_solved = True
except:
    scamp_solved = False
try:
    an_solved = str(mosaic_header["A_ORDER"])
    if an_solved > 0:
        an_solved = True
except:
    an_solved = False
mosaic_hdulist.close()
if scamp_solved or an_solved:
    system(sethead_bin + " h_long_" + obs_string + "_coadd.fits WCSFIT=True")
    system(sethead_bin + " h_long_" + obs_string + "_coadd.weight.fits WCSFIT=True")
    if do_short:
        system(sethead_bin + " h_short_" + obs_string + "_coadd.fits WCSFIT=True")
        system(sethead_bin + " h_short_" + obs_string + "_coadd.weight.fits WCSFIT=True")
else:
    system(sethead_bin + " h_long_" + obs_string + "_coadd.fits WCSFIT=False")
    try:
        system(python_bin + " autoastrometry.py h_long_" + obs_string + "_coadd.fits -px 1.00 -inv -c ub2")
        h_long_hdulist = pyfits.open("h_long_" + obs_string + "_coadd.fits", "readonly")
        h_long_weights_hdulist = pyfits.open("h_long_" + obs_string + 
            "_coadd.weight.fits", "update")
        if do_short:
            h_short_hdulist = pyfits.open("h_short_" + obs_string + "_coadd.fits", 
                "update")
            h_short_weights_hdulist = pyfits.open("h_short_" + obs_string + 
                "_coadd.weight.fits", "update")
        tmp0 = h_long_hdulist[0].header
        del tmp0["SIMPLE"]
        del tmp0["BITPIX"]
        del tmp0["NAXIS"]
        h_long_hdulist.close()
        for c in tmp0.ascardlist():
            h_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
            if do_short:
                h_short_hdulist[0].header.update(c.key,c.value,c.comment)
                h_short_weights_hdulist[0].header.update(c.key,c.value,c.comment)
        h_long_weights_hdulist.verify("silentfix")
        h_long_weights_hdulist.close(output_verify='warn')
        if do_short:
            h_short_hdulist.verify("silentfix")
            h_short_hdulist.close(output_verify='warn')
            h_short_weights_hdulist.verify("silentfix")
            h_short_weights_hdulist.close(output_verify='warn')
    except:
        print "Failed H WCS fit."

mosaic_hdulist = pyfits.open("k_long_" + obs_string + "_coadd.fits")
mosaic_header = mosaic_hdulist[0].header
try:
    scamp_solved = float(mosaic_header["PV1_1"])
    if scamp_solved > 0:
        scamp_solved = True
except:
    scamp_solved = False
try:
    an_solved = str(mosaic_header["A_ORDER"])
    if an_solved > 0:
        an_solved = True
except:
    an_solved = False
mosaic_hdulist.close()
if scamp_solved or an_solved:
    system(sethead_bin + " k_long_" + obs_string + "_coadd.fits WCSFIT=True")
    system(sethead_bin + " k_long_" + obs_string + "_coadd.weight.fits WCSFIT=True")
    if do_short:
        system(sethead_bin + " k_short_" + obs_string + "_coadd.fits WCSFIT=True")
        system(sethead_bin + " k_short_" + obs_string + "_coadd.weight.fits WCSFIT=True")
else:
    system(sethead_bin + " k_long_" + obs_string + "_coadd.fits WCSFIT=False")
    try:
        system(python_bin + " autoastrometry.py k_long_" + obs_string + "_coadd.fits -px 1.00 -inv -c ub2")
        k_long_hdulist = pyfits.open("k_long_" + obs_string + "_coadd.fits", "readonly")
        k_long_weights_hdulist = pyfits.open("k_long_" + obs_string + 
            "_coadd.weight.fits", "update")
        if do_short:
            k_short_hdulist = pyfits.open("k_short_" + obs_string + "_coadd.fits", 
                "update")
            k_short_weights_hdulist = pyfits.open("k_short_" + obs_string + 
                "_coadd.weight.fits", "update")
        tmp0 = k_long_hdulist[0].header
        del tmp0["SIMPLE"]
        del tmp0["BITPIX"]
        del tmp0["NAXIS"]
        k_long_hdulist.close()
        for c in tmp0.ascardlist():
            k_long_weights_hdulist[0].header.update(c.key,c.value,c.comment)
            if do_short:
                k_short_hdulist[0].header.update(c.key,c.value,c.comment)
                k_short_weights_hdulist[0].header.update(c.key,c.value,c.comment)
        k_long_weights_hdulist.verify("silentfix")
        k_long_weights_hdulist.close(output_verify='warn')
        if do_short:
            k_short_hdulist.verify("silentfix")
            k_short_hdulist.close(output_verify='warn')
            k_short_weights_hdulist.verify("silentfix")
            k_short_weights_hdulist.close(output_verify='warn')
    except:
        print "Failed K WCS fit."
# Set the NSWARP header keyword.
nswarp_number = num_dither_positions * 3
system(sethead_bin + " j_long_" + obs_string + "_coadd.fits NSWARP=" + 
    str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
system(sethead_bin + " h_long_" + obs_string + "_coadd.fits NSWARP=" + 
    str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
system(sethead_bin + " k_long_" + obs_string + "_coadd.fits NSWARP=" + 
    str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
system(sethead_bin + " j_long_" + obs_string + "_coadd.weight.fits NSWARP=" + 
    str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
system(sethead_bin + " h_long_" + obs_string + "_coadd.weight.fits NSWARP=" + 
    str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
system(sethead_bin + " k_long_" + obs_string + "_coadd.weight.fits NSWARP=" + 
    str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
if do_short:
    system(sethead_bin + " j_short_" + obs_string + "_coadd.fits NSWARP=" + 
        str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
    system(sethead_bin + " h_short_" + obs_string + "_coadd.fits NSWARP=" + 
        str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
    system(sethead_bin + " k_short_" + obs_string + "_coadd.fits NSWARP=" + 
        str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
    system(sethead_bin + " j_short_" + obs_string + "_coadd.weight.fits NSWARP=" + 
        str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
    system(sethead_bin + " h_short_" + obs_string + "_coadd.weight.fits NSWARP=" + 
        str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
    system(sethead_bin + " k_short_" + obs_string + "_coadd.weight.fits NSWARP=" + 
        str(nswarp_number) + ''' / "num files in used to make this mosaic"''')
# Set the RA and DEC header keywords in the final mosaics.
mosaic_hdulist = pyfits.open("j_long_" + obs_string + "_coadd.fits")
mosaic_header = mosaic_hdulist[0].header
crval1 = float(mosaic_header["CRVAL1"])
crval2 = float(mosaic_header["CRVAL2"])
crpix1 = float(mosaic_header["CRPIX1"])
crpix2 = float(mosaic_header["CRPIX2"])
mosaic_hdulist.close()
mosaic_dec = crval2 + (256.0 - crpix2)/3600
mosaic_ra = crval1 + (crpix1 - 256.0)/3600/cos(mosaic_dec*0.0174532925)
mosaic_dec = str(ephem.degrees(mosaic_dec * 0.0174532925))
mosaic_ra = str(ephem.hours(mosaic_ra * 0.0174532925))
system(sethead_bin + " j_long_" + obs_string + "_coadd.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " j_long_" + obs_string + "_coadd.weight.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " j_long_" + obs_string + "_coadd.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
system(sethead_bin + " j_long_" + obs_string + "_coadd.weight.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
if do_short:
    system(sethead_bin + " j_short_" + obs_string + "_coadd.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " j_short_" + obs_string + "_coadd.weight.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " j_short_" + obs_string + "_coadd.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
    system(sethead_bin + " j_short_" + obs_string + "_coadd.weight.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')

mosaic_hdulist = pyfits.open("h_long_" + obs_string + "_coadd.fits")
mosaic_header = mosaic_hdulist[0].header
crval1 = float(mosaic_header["CRVAL1"])
crval2 = float(mosaic_header["CRVAL2"])
crpix1 = float(mosaic_header["CRPIX1"])
crpix2 = float(mosaic_header["CRPIX2"])
mosaic_hdulist.close()
mosaic_dec = crval2 + (256.0 - crpix2)/3600
mosaic_ra = crval1 + (crpix1 - 256.0)/3600/cos(mosaic_dec*0.0174532925)
mosaic_dec = str(ephem.degrees(mosaic_dec * 0.0174532925))
mosaic_ra = str(ephem.hours(mosaic_ra * 0.0174532925))
system(sethead_bin + " h_long_" + obs_string + "_coadd.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " h_long_" + obs_string + "_coadd.weight.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " h_long_" + obs_string + "_coadd.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
system(sethead_bin + " h_long_" + obs_string + "_coadd.weight.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
if do_short:
    system(sethead_bin + " h_short_" + obs_string + "_coadd.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " h_short_" + obs_string + "_coadd.weight.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " h_short_" + obs_string + "_coadd.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
    system(sethead_bin + " h_short_" + obs_string + "_coadd.weight.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
        
mosaic_hdulist = pyfits.open("k_long_" + obs_string + "_coadd.fits")
mosaic_header = mosaic_hdulist[0].header
crval1 = float(mosaic_header["CRVAL1"])
crval2 = float(mosaic_header["CRVAL2"])
crpix1 = float(mosaic_header["CRPIX1"])
crpix2 = float(mosaic_header["CRPIX2"])
mosaic_hdulist.close()
mosaic_dec = crval2 + (256.0 - crpix2)/3600
mosaic_ra = crval1 + (crpix1 - 256.0)/3600/cos(mosaic_dec*0.0174532925)
mosaic_dec = str(ephem.degrees(mosaic_dec * 0.0174532925))
mosaic_ra = str(ephem.hours(mosaic_ra * 0.0174532925))
system(sethead_bin + " k_long_" + obs_string + "_coadd.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " k_long_" + obs_string + "_coadd.weight.fits RA=" + 
    str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
system(sethead_bin + " k_long_" + obs_string + "_coadd.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
system(sethead_bin + " k_long_" + obs_string + "_coadd.weight.fits DEC=" + 
    str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
if do_short:
    system(sethead_bin + " k_short_" + obs_string + "_coadd.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " k_short_" + obs_string + "_coadd.weight.fits RA=" + 
        str(mosaic_ra) + ''' / "right ascension coordinate of mosaic center"''')
    system(sethead_bin + " k_short_" + obs_string + "_coadd.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')
    system(sethead_bin + " k_short_" + obs_string + "_coadd.weight.fits DEC=" + 
        str(mosaic_dec) + ''' / "declination coordinate of mosaic center"''')

# Set the SATNOLN header keyword for the mosaics.
system(sethead_bin + " j_long_" + obs_string + "_coadd.fits SATNONLN=" + 
    '''5382000 / "Saturation due to CCD non-linearity"''')
system(sethead_bin + " h_long_" + obs_string + "_coadd.fits SATNONLN=" + 
    '''5382000 / "Saturation due to CCD non-linearity"''')
system(sethead_bin + " k_long_" + obs_string + "_coadd.fits SATNONLN=" + 
    '''5382000 / "Saturation due to CCD non-linearity"''')
system(sethead_bin + " j_long_" + obs_string + "_coadd.weight.fits SATNONLN=" + 
    '''5382000 / "Saturation due to CCD non-linearity"''')
system(sethead_bin + " h_long_" + obs_string + "_coadd.weight.fits SATNONLN=" + 
    '''5382000 / "Saturation due to CCD non-linearity"''')
system(sethead_bin + " k_long_" + obs_string + "_coadd.weight.fits SATNONLN=" + 
    '''5382000 / "Saturation due to CCD non-linearity"''')
if do_short:
    system(sethead_bin + " j_short_" + obs_string + "_coadd.fits SATNONLN=" + 
        '''5382000 / "Saturation due to CCD non-linearity"''')
    system(sethead_bin + " h_short_" + obs_string + "_coadd.fits SATNONLN=" + 
        '''5382000 / "Saturation due to CCD non-linearity"''')
    system(sethead_bin + " k_short_" + obs_string + "_coadd.fits SATNONLN=" + 
        '''5382000 / "Saturation due to CCD non-linearity"''')
    system(sethead_bin + " j_short_" + obs_string + "_coadd.weight.fits SATNONLN=" + 
        '''5382000 / "Saturation due to CCD non-linearity"''')
    system(sethead_bin + " h_short_" + obs_string + "_coadd.weight.fits SATNONLN=" + 
        '''5382000 / "Saturation due to CCD non-linearity"''')
    system(sethead_bin + " k_short_" + obs_string + "_coadd.weight.fits SATNONLN=" + 
        '''5382000 / "Saturation due to CCD non-linearity"''')
        
print "Finished wcs fitting coadd images, %f seconds required." % (time() - t1)
# Clean up the directory of all the intermediate files.
system("rm -rf " + obs_string + "_raw")
file_date = strt_cpu.split(" ")[0]
system("mkdir " + obs_string + "_" + file_date + "-reduction_output")
if do_full:
    system("mv -f weightmaps " + obs_string + "_weightmaps")
    system("mv -f triplestacks " + obs_string + "_triplestacks")
    system("mv -f triplestackweights " + obs_string + "_triplestackweights")
    system("mv -f reduced " + obs_string + "_reduced")
system("mkdir ./" + obs_string + "_mosaics")
system("mv -f *" + obs_string + "_coadd*fits ./" + obs_string + "_mosaics")
system("mv *" + obs_string + "_* " + obs_string + "_" + file_date + "-reduction_output")
system("mv ./*alignments*.txt " + obs_string + "_" + file_date + "-reduction_output")
system("mkdir " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + 
    "_ancillary")
system("mv *_triplestack*.txt " + obs_string + "_" + file_date + "-reduction_output/" + 
    obs_string + "_ancillary")
system("mv cloud_rejection.txt " + obs_string + "_" + file_date + "-reduction_output/" + 
    obs_string + "_ancillary")
system("mv sethead_commands.txt " + obs_string + "_" + file_date + "-reduction_output/" + 
    obs_string + "_ancillary")
system("mkdir " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + 
    "_alignments")
system("mv " + obs_string + "_" + file_date + "-reduction_output/*.txt " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + "_alignments")
if DEBUG:
    system("mkdir " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + 
        "_debug")
    system("mv gaussianfiltered " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + "_debug")
    system("mv dynamic_bpm " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + "_debug") 
    system("mv finalweights " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + "_debug")
    system("mv sigmagaussian " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + "_debug")
    system("mv objectimages " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + "_debug")
    system("mv masks " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + 
        "_debug")
    system("mv skarks " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + 
        "_debug")
    system("mv midreduced_triplestackweights " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + "_debug")
    system("mv midreduced_triplestacks " + obs_string + "_" + file_date + "-reduction_output/" + obs_string + "_debug")
system("scp -P 10222 -qr " + obs_string + "_" + file_date + "-reduction_output " + dest_directory)
system("rm -rf /home/ptelreducer/scratch/" + obs_string + "-reduction_directory")
# End program execution timing.
end_time = time()
total_time = end_time - start_time
print "Program finished, execution time %f seconds." % total_time
if astromery_dot_net_worked == False:
    print "Warning: astronometry.net failed for the triplestack sourceimage."
    print "    This means WCS of triplestacks and reduced images is inaccurate."
