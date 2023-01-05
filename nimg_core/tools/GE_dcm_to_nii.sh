#!/bin/bash

###########################################################################
# Copyright (c) 2017 Vincent A. Magnotta & Joel Bruss
#   With a strong assist from James Kent!
# All rights reserved.
#
# Redistribution and use in source and binary forms are permitted
# provided that the above copyright notice and this paragraph are
# duplicated in all such forms and that any documentation,
# advertising materials, and other materials related to such
# distribution and use acknowledge that the software was developed
# by the University of Iowa MR Research Facility (MRRF). The name of the
# MRRF may not be used to endorse or promote products derived
# from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
###########################################################################

###########################################################################
#
# Script to Convert fMRI images from a GE scanner to 4D image format (NIFTI
# or AFNI). This script has been tested on images from GE scanners running
# DV23 and DV25 software versions.
#
# This script needs AFNI, FSL, FreeSurfer and gdcm to run correctly
#  AFNI:  https://afni.nimh.nih.gov/afni/
#  FSL:  https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL
#  FreeSurfer:  https://surfer.nmr.mgh.harvard.edu/
#  gdcm:  http://gdcm.sourceforge.net/html/gdcmdump.html
#
# Darwin (Mac) may also require Homebrew (for gdcm):
#  Howebrew: https://brew.sh/
#
#  This script only allow for EPI data to be reconstructed in a sequential
#  or interleaved pattern, using either seq+z or alt+z via AFNI's 'to3d'.
#  Scan collection, with an oppopsite reconstruction (alt/seq-z) is not
#  handled.  Further, GE only allos acquisition starting at slice "0", so
#  one could never have alt+z2,
###########################################################################

scriptdir=`dirname $(perl -e 'use Cwd "abs_path";print abs_path(shift)' $0)`
dcmVerbose=0
reorientFlag=0
imageExt=".nii.gz"
imageType="NIfTI"

function printCommandLine {
  echo ""
  echo "Usage: GE_dcm_to_nii.sh -i DICOM_directory -o outputDirectory -b outputBaseName -v -r -a -B"
  echo ""
  echo "  e.g.:"
  echo "       GE_dcm_to_nii.sh -i /path/to/DICOM -o /output/directory -b fmri_Visual -v -r -a -B"
  echo ""
  echo "   where:"
  echo "   -i Directory where DICOM data is located"
  echo "   -o Output directory"
  echo "       *If this is not specified, this defaults to input DICOM directory"
  echo "   -b Output basename (e.g. fmri_Visual)"
  echo "       *If this is not specified, this will be pulled from 'ID Series Description' in the DICOM header"
  echo "   -v Verbose output"
  echo "       *Additional file information will be logged to Output directory"
  echo "         *This only applies to EPI images"
  echo "           *This option is off by default"
  echo "   -r Reorient file to RPI (FSL's MNI preferred orientation)"
  echo "       *This option is off by default"
  echo "         *This will be skipped for localizers"
  echo "   -a AFNI BRIKs"
  echo "       *Output AFNI BRIKs instead of NIFTI formatted images"
  echo "       *NIFTI image is the default output type"
  echo "   -B Output a BIDS compatible .json file"
  echo ""
  exit 1
}


# Parse Command line arguments
while getopts “hi:o:b:Bvra” OPTION
do
  case $OPTION in
    h)
      printCommandLine
      ;;
    i)
      dcmDir=$OPTARG
      ;;
    o)
      outDir=$OPTARG
      ;;
    b)
      outBase=$OPTARG
      ;;
    v)
      dcmVerbose=1
      ;;
    r)
      reorientFlag=1
      ;;
    a)
      imageExt=""
      imageType="AFNI"
      ;;
    B)
      BIDS=true
      ;;
    ?)
      printCommandLine
      ;;
     esac
done

##################
#function: native_orient
##################
#purpose: returns filename without it's directory path or file extensions
##################
#dependencies: fsl
##################
#Used in: main
##################
function native_orient() {
    local infile=$1 &&\
    [ ! -z  "${infile}" ] ||\
    ( printf '%s\n' "${FUNCNAME[0]}, input not defined" && return 1 )

    #Determine qform-orientation to properly reorient file to RPI (MNI) orientation
  xorient=`fslhd ${infile} | grep "^qform_xorient" | awk '{print $2}' | cut -c1`
  yorient=`fslhd ${infile} | grep "^qform_yorient" | awk '{print $2}' | cut -c1`
  zorient=`fslhd ${infile} | grep "^qform_zorient" | awk '{print $2}' | cut -c1`

  native_orient=${xorient}${yorient}${zorient}
  echo ""
  echo "native orientation = ${native_orient}"
}

##################
#function: RPI_orient
##################
#purpose: returns filename without it's directory path or file extensions
##################
#dependencies: fsl
##################
#Used in: main
##################
function RPI_orient() {
    local infile=$1 &&\
    [ ! -z  "${infile}" ] ||\
    ( printf '%s\n' "${FUNCNAME[0]}, input not defined" && return 1 )

    #Determine qform-orientation to properly reorient file to RPI (MNI) orientation
  #xorient=`fslhd ${infile} | grep "^qform_xorient" | awk '{print $2}' | cut -c1`
  #yorient=`fslhd ${infile} | grep "^qform_yorient" | awk '{print $2}' | cut -c1`
  #zorient=`fslhd ${infile} | grep "^qform_zorient" | awk '{print $2}' | cut -c1`

  #native_orient=${xorient}${yorient}${zorient}

  #echo "native orientation = ${native_orient}"

  if [ "${native_orient}" != "RPI" ]; then

    case ${native_orient} in

    #L PA IS
    LPI)
      flipFlag="-x y z"
      ;;
    LPS)
      flipFlag="-x y -z"
      ;;
    LAI)
      flipFlag="-x -y z"
      ;;
    LAS)
      flipFlag="-x -y -z"
      ;;

    #R PA IS
    RPS)
      flipFlag="x y -z"
      ;;
    RAI)
      flipFlag="x -y z"
      ;;
    RAS)
      flipFlag="x -y -z"
      ;;

    #L IS PA
    LIP)
      flipFlag="-x z y"
      ;;
    LIA)
      flipFlag="-x -z y"
      ;;
    LSP)
      flipFlag="-x z -y"
      ;;
    LSA)
      flipFlag="-x -z -y"
      ;;

    #R IS PA
    RIP)
      flipFlag="x z y"
      ;;
    RIA)
      flipFlag="x -z y"
      ;;
    RSP)
      flipFlag="x z -y"
      ;;
    RSA)
      flipFlag="x -z -y"
      ;;

    #P IS LR
    PIL)
      flipFlag="-z x y"
      ;;
    PIR)
      flipFlag="z x y"
      ;;
    PSL)
      flipFlag="-z x -y"
      ;;
    PSR)
      flipFlag="z x -y"
      ;;

    #A IS LR
    AIL)
      flipFlag="-z -x y"
      ;;
    AIR)
      flipFlag="z -x y"
      ;;
    ASL)
      flipFlag="-z -x -y"
      ;;
    ASR)
      flipFlag="z -x -y"
      ;;

    #P LR IS
    PLI)
      flipFlag="-y x z"
      ;;
    PLS)
      flipFlag="-y x -z"
      ;;
    PRI)
      flipFlag="y x z"
      ;;
    PRS)
      flipFlag="y x -z"
      ;;

    #A LR IS
    ALI)
      flipFlag="-y -x z"
      ;;
    ALS)
      flipFlag="-y -x -z"
      ;;
    ARI)
      flipFlag="y -x z"
      ;;
    ARS)
      flipFlag="y -x -z"
      ;;

    #I LR PA
    ILP)
      flipFlag="-y z x"
      ;;
    ILA)
      flipFlag="-y -z x"
      ;;
    IRP)
      flipFlag="y z x"
      ;;
    IRA)
      flipFlag="y -z x"
      ;;

    #S LR PA
    SLP)
      flipFlag="-y z -x"
      ;;
    SLA)
      flipFlag="-y -z -x"
      ;;
    SRP)
      flipFlag="y z -x"
      ;;
    SRA)
      flipFlag="y -z -x"
      ;;

    #I PA LR
    IPL)
      flipFlag="-z y x"
      ;;
    IPR)
      flipFlag="z y x"
      ;;
    IAL)
      flipFlag="-z -y x"
      ;;
    IAR)
      flipFlag="z -y x"
      ;;

    #S PA LR
    SPL)
      flipFlag="-z y -x"
      ;;
    SPR)
      flipFlag="z y -x"
      ;;
    SAL)
      flipFlag="-z -y -x"
      ;;
    SAR)
      flipFlag="z -y -x"
      ;;
    esac

    echo "flipping by ${flipFlag}"

    #Reorienting image and checking for warning messages
    warnFlag=`fslswapdim ${infile} ${flipFlag} ${infile%.nii.gz}_RPI.nii.gz`
    warnFlagCut=`echo ${warnFlag} | awk -F":" '{print $1}'`

    #Reorienting the file may require swapping out the flag orientation to match the .img block
    if [[ $warnFlagCut == "WARNING" ]]; then
      fslorient -swaporient ${infile%.nii.gz}_RPI.nii.gz
    fi

    #Remove the intermediate file
    rm $infile

  else
    #Data already in RPI orientation
    echo "No need to reorient.  Dataset already in RPI orientation."

    if [ ! -e ${infile%.nii.gz}_RPI.nii.gz ]; then
      mv ${infile} ${infile%.nii.gz}_RPI.nii.gz
    fi

  fi
}

##################
#function: cpu_thread_num
##################
#purpose: finds the number of processors (logical) on your machine
##################
#dependencies: none
##################
#Used in: main
##################
function cpu_thread_num(){
  local arch=$(uname)
  if [ "${arch}" == "Darwin" ];then
    local cpu_threads=$(sysctl -n hw.ncpu)
  elif [ "${arch}" == "Linux" ];then
    local cpu_threads=$(nproc)
  else
    printf "WARNING: I don't know what OS you are running, assuming Linux" 1>&2
    local cpu_threads=$(nproc)
  fi
  printf "${cpu_threads}"
}


#Check for AFNI
  #Used for conversion to NIfTI, reorienting to RPI (if chosen), clipping first TR (multiband)
dcmHdrProg=`which dicom_hdr`
if [[ $dcmHdrProg == "" ]]; then
  echo "Error:  Unable to find the AFNI DICOM utility (dicom_hdr). Update your path and rerun the command."
  exit 1
fi

#Check for gdcmdump
  #Verified that this can be installed via homebrew (http://brew.sh/) for MacOS
gdcmProg=`which gdcmdump`
if [[ $gdcmProg == "" ]]; then
  echo "Error:  Unable to find gdcmdump. Update your path and rerun the command."
  echo "  Note:  If using MacOS, you *should* be able to install this with homebrew."
  exit 1
fi

#Check for fslswapdim, ONLY used to reorient NIfTI output
if [[ $imageType == "NIFTI" ]]; then
  if [[ $reorientFlag -eq 1 ]]; then
    fslProg=`which fslswapdim`
    if [[ $fslProg == "" ]]; then
      echo "Error:  Unable to find the FSL utility (fslswapdim). Update your path and rerun the command."
      exit 1
    fi
  fi
fi

#Check for mri_convert (FreeSurfer), used to convert anything other than EPI, loclaizer or calibration scans
fsProg=`which mri_convert`
if [[ $fsProg == "" ]]; then
  echo "Error:  Unable to find the FreeSurfer DICOM utility (mri_convert). Update your path and rerun the command."
  exit 1
fi

#If options aren't set at input, default to the following
if [[ $dcmDir == "" ]]; then
  echo "Error:  No input directory chosen.  Please point to a dirctory where DICOM data is located with the '-i' option"
  exit 1
fi

if [[ $outDir == "" ]]; then
  outDir=${dcmDir}
else
  #If $outDir doesn't exist, make it
  if [[ ! -d $outDir ]]; then
    mkdir -p $outDir
  fi
fi


#pushd `pwd` "$@" > /dev/null
pushd $PWD > /dev/null



#######################################################################
# Begin pulling in information about the DICOM series
#######################################################################

echo ""
echo "Gathering DICOM parameters"
echo ""

pushd ${dcmDir} > /dev/null

#Source a DICOM file to strip header information from
dcmList=($(find . -maxdepth 1 -name '*.dcm' | sed -e "s|^./||g" ))
dcmPic=${dcmList[1]}

if [[ `echo $dcmPic` == "" ]]; then
  echo "Error:  No DICOM is present in the input directory set with '-i.'  Please check the PATH and try again."
  exit 1
fi

#Determining information about scan (number of volumes (TRs) and number of slices per volume, TR and TE, etc.)
softRev=`dicom_hdr $dcmPic | grep "ACQ Software Version" | awk -F ":" '{print $2}' | awk -F "." '{print $1}'`
scanMan=`dicom_hdr $dcmPic | grep "ID Manufacturer/" | awk -F"//" '{print $3}' | tr -s ' ' | tr ' ' '_'`
scanModel=`dicom_hdr $dcmPic | grep "ID Manufacturer Model Name" | awk -F"//" '{print $3}' | tr -s ' ' | tr ' ' '_'`
magSize=`dicom_hdr $dcmPic | grep "ACQ Magnetic Field Strength" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
dcmTR=`dicom_hdr $dcmPic | grep "0018 0080" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
dcmTE=`dicom_hdr $dcmPic | grep "0018 0081" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
sliceThickness=`dicom_hdr $dcmPic | grep "ACQ Slice Thickness" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
flipAng=`dicom_hdr $dcmPic | grep "ACQ Flip Angle" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
matX=`gdcmdump -P -i $dcmPic | grep "MATRIXX" | awk -F'"' '{print $2}'`
matY=`gdcmdump -P -i $dcmPic | grep "MATRIXY" | awk -F'"' '{print $2}'`
matrixSize="${matX}x${matY}"
#Number of prescribed frames
dcmVolTot=`dicom_hdr $dcmPic | grep "0020 0105" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
dcmTot=`dicom_hdr $dcmPic | grep "0020 1002" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
dcmRelSeries=`dicom_hdr $dcmPic | grep "0020 0011" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
#important for BIDs determining format
dcmSeries=`dicom_hdr $dcmPic | grep "0008 103e" | awk -F"//" '{print $3}' | tr -d '[[:space:]]' | sed 's|'/'|'_'|g' | sed 's|':'|'_'|g'`
dcmSubSeries=`dicom_hdr $dcmPic | grep "0019 109c" | awk -F"//" '{print $3}' | tr -d '[[:space:]]' | sed 's|'/'|'_'|g' | sed 's|':'|'_'|g'`
  #Now have conditions where "epiRT" could mean full functional run or 1000+ series "mean" versions
tmp_dcmVols=`dicom_hdr $dcmPic | grep "0043 1079" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`


  #Need to adjust for multiband/EPI
    #Remove 5 volumes (calibration data)  (One more volume, post conversion, will be stripped off later) (v1.0 multiband)
#if [[ $dcmSubSeries == "multiband_mux_epi" ]]; then
  #dcmNumVols=`echo $dcmVolTot | awk '{print ($1-5)}'`
  #dcmNumSlices=`dicom_hdr $dcmPic | grep "0021 104f" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
if [[ $dcmSubSeries == "multiband_mux_epi" ]]; then
  if [[ $softRev == "DV26" ]]; then
    dcmNumBands=`gdcmdump -P -i $dcmPic | grep "MBACCEL" | awk -F"\"" '{print $2}'`
    dcmNumVols=$dcmVolTot
    dcmNumSlices=`dicom_hdr $dcmPic | grep "0021 104f" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
  else
    dcmNumBands=`gdcmdump -P -i $dcmPic | grep "USERCV22" | awk -F"\"" '{print $2}'`
    dcmNumVols=$dcmVolTot
    dcmNumSlices=`dicom_hdr $dcmPic | grep "0021 104f" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
  fi
elif [[ $dcmSubSeries == "multiband_epi" ]]; then
#if [[ $dcmSubSeries == "multiband_mux_epi" || $dcmSubSeries == "multiband_epi" ]]; then
  ##As of 8/2017, the formula for "true" number of volumes is as follows:
    #${dcmVolTot}-(2*# of bands)+1 avg. calibration frame
  ##11/2017 update.  New version of multiband that re-uses "multband_mux_epi" but is defined the same as "multiband_epi"
    #In conversations with Vince Magnotta, we're under the strong assumption that no one is using the really, really old multiband sequence.  That conversion method is commented out.
  #Number of bands (currently a user-defined field (22), but will eventually have a dedicated GE DICOM field)
  dcmNumBands=`gdcmdump -P -i $dcmPic | grep "USERCV22" | awk -F"\"" '{print $2}'`
  dcmNumVols=`echo $dcmVolTot $dcmNumBands | awk '{print ($1-($2*2)+1)}'`
  dcmNumSlices=`dicom_hdr $dcmPic | grep "0021 104f" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
elif [[ $dcmSubSeries == "epiRT" || $dcmSubSeries == "ia_stable_epiRTTopUp" ]]; then
  dcmNumSlices=`dicom_hdr $dcmPic | grep "0043 1079" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
  if [[ $magSize -eq 7 ]]; then
    dcmNumVols=$dcmVolTot
  else
    dcmNumVols=`dicom_hdr $dcmPic | grep "2001 1041" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
  fi
elif [[ $dcmSubSeries == "ssfse" ]]; then
  dcmNumSlices=`dicom_hdr $dcmPic | grep "0021 104f" | awk -F"//" '{print $3}' | tr -d '[[:space:]]' | awk '{print ($1/3)}'`
  tmpdcmNumVols=`dicom_hdr $dcmPic | grep "0020 1002" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
  dcmNumVols=`echo $tmpdcmNumVols $dcmNumSlices | awk '{print ($1/$2)}'`
else
  dcmNumVols=1
  dcmNumSlices=`dicom_hdr $dcmPic | grep "0020 1002" | awk -F"//" '{print $3}' | tr -d '[[:space:]]'`
fi


##Determine Phase encoding direction (for EPI images), dwell Time (s)
if [[ $dcmSubSeries == "multiband_mux_epi" || $dcmSubSeries == "epiRT" || $dcmSubSeries == "multiband_epi" || $dcmSubSeries == "ia_stable_epiRTTopUp" ]]; then

  ### Slice Orientation #########

  #Pull out the Primary/Secondary orientation information from DICOM (will use to determine if coronal/axial/sagittal acquisition)
  tmpdcmDir=`dicom_hdr $dcmPic | grep "0020 0037" | awk -F"//" '{print $NF}' | tr -d '[[:space:]]'`

  #Convert the "\", possibly obliqued primary/secondary vector information into one long, six-digit string of absolute values, cardinal orientations
    #e.g. go from "1\-0\0\-0\1\0" or "0.985\-0\0\-0.324\0.89\0" to "100010"
    #Note:  This *should* convert andy oblique dataset to a cardinal set, but I haven't been able to test on many cases.  There is a chance that a highly angled, e.g. axial will show up as a coronal due to rounding. (JB)
  primSeconVector=`paste -d '\t' -s <(echo $tmpdcmDir | awk -F'\' '{for(i=1;i<=NF;i++){printf ("%3.0f\n",$i)}}' | sed 's|-||g') | sed 's/[[:blank:]]//g'`

  #Use the six-digit Primary/Secondary vector code to determine X/Y/Z orientation axes
  if [[ $primSeconVector -eq 100010 ]]; then
    primSeconDir=XY
  elif [[ $primSeconVector -eq 100001 ]]; then
    primSeconDir=XZ
  elif [[ $primSeconVector -eq 010001 ]]; then
    primSeconDir=YZ
  else
    echo "${primSeconVector} is an unkown Primary/Secondary ordering.  Unable to process EPI data."
    exit 1
  fi

  #Combine the Primary and Secondary Axes, determine the layout to be Coronal, Axial or Sagittal
  #scanVectors=${primDir}${seconDir}

  #Use the Primary and Secondary Vectors to determine scan type
  case "${primSeconDir}" in
    XY)
      scanType=Axial
      ;;
    XZ)
      scanType=Coronal
      ;;
    YZ)
      scanType=Sagittal
      ;;
    *)
      scanType=Unknown
      echo "Unknown acquisition ${scanVectors} from a '${primDir},${seconDir}' orientation."
      exit 1
  esac

  ###############################


  ### Dwell Time ################

  #The "true" dwell Time can be calculated from the dwellTime reported in the DICOM header divided by the Acceleration Factor
    #Acceleration Factor = 1/Asset Factor (from DICOM) -OR- pulled directly from the compressed DICOM field 0025,1b (PHASEACCEL):
      #https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1701&L=FSL&F=&S=&P=16997
      #https://groups.google.com/forum/#!topic/comp.protocols.dicom/mxnCkv8A-i4

  #Pull the dwellTime, uncorrected for Acceleration Factor, from the DICOM header, convert from micro-seconds to seconds
  tmpDwell=`dicom_hdr $dcmPic | grep "0043 102c" | awk -F"//" '{print $NF}' | tr -d '[[:space:]]' | awk '{print $1/1000000}'`
  echoSpacing=`echo $tmpDwell | awk '{printf "%f", ($1*1000)}'`
  #Pull out the Asset Factor
    #Acceleration Factor = 1/Asset Factor
      #-OR- `gdcmdump -P -i $dcmPic | grep "PHASEACCEL" | awk -F'"' '{print $2}'`
  assetFactor=`dicom_hdr $dcmPic | grep "0043 1083" | awk -F"//" '{print $NF}' | tr -d '[[:space:]]'`
  accelFactor=`echo $assetFactor | awk '{print 1/$1}'`
  #Correct the temp dwellTime for the Acceleration Factor, to get dwellTime in S
  dwellTime=`echo $tmpDwell $accelFactor | awk '{print $1/$2}'`

  ###############################


  ### Phase Encoding direction ##
    #http://fmri.ucsd.edu/Howto/3T/fmapproc.html

  #Look for ROW/COL phase encoding direction
  dcmOrient=`dicom_hdr $dcmPic | grep "0018 1312" | awk -F"//" '{print $NF}' | tr -d '[[:space:]]'`

  if [[ $scanType == "Axial" ]]; then
    if [[ $dcmOrient == "COL" ]]; then
      peDir=AP
    elif [[ $dcmOrient == "ROW" ]]; then
      peDir=RL
    else
      peDir=Unknown
    fi
  elif [[ $scanType == "Sagittal" ]]; then
    if [[ $dcmOrient == "COL" ]]; then
      peDir=SI
    elif [[ $dcmOrient == "ROW" ]]; then
      peDir=AP
    else
      peDir=Unknown
    fi
  elif [[ $scanType == "Coronal" ]]; then
    if [[ $dcmOrient == "COL" ]]; then
      peDir=SI
    elif [[ $dcmOrient == "ROW" ]]; then
      peDir=RL
    else
      peDir=Unknown
    fi
  else
    peDir=Unknown
  fi

  ###############################


  ### View Order ################

    #VO=0=Top-Down
    #VO=1=Bottom-Up

  tmpviewOrder=`gdcmdump -P $dcmPic | grep "VIEWORDER" | awk '{print $2}' | sed 's/\"//g'`
  if [[ $tmpviewOrder -eq 0 ]]; then
    viewOrder=TD
  else
    viewOrder=BU
  fi

  ###############################


  ### unwarpdir #################

    #Put Slice Orientation, Phase Encoding Direction and View order together to determine "unwarpdir" (to be used with, along with dwellTime, "epi_reg")
    #http://fmri.ucsd.edu/Howto/3T/fmapproc.html

    #epi_reg, although demands the "sign" of the unwarpdir first (e.g. -y), it actually sets it opposite, internally

  if [[ $scanType == "Axial" ]]; then
    if [[ $peDir == "AP" ]]; then
      if [[ $viewOrder == "BU" ]]; then
        unWarpDir="y"
      elif [[ $viewOrder == "TD" ]]; then
        unWarpDir="-y"
      else
        unWarpDir=Unknown
      fi
    elif [[ $peDir == "RL" ]]; then
      if [[ $viewOrder == "BU" ]]; then
        unWarpDir="-x"
      elif [[ $viewOrder == "TD" ]]; then
        unWarpDir="x"
      else
        unWarpDir=Unknown
      fi
    else
      unWarpDir=Unknown
    fi

  elif [[ $scanType == "Sagittal" ]]; then
    if [[ $peDir == "SI" ]]; then
      if [[ $viewOrder == "BU" ]]; then
        unWarpDir="-y"
      elif [[ $viewOrder == "TD" ]]; then
        unWarpDir="y"
      else
        unWarpDir=Unknown
      fi
    elif [[ $peDir == "AP" ]]; then
      if [[ $viewOrder == "BU" ]]; then
        unWarpDir="x"
      elif [[ $viewOrder == "TD" ]]; then
        unWarpDir="-x"
      else
        unWarpDir=Unknown
      fi
    else
      unWarpDir=Unknown
    fi

  elif [[ $scanType == "Coronal" ]]; then
    if [[ $peDir == "SI" ]]; then
      if [[ $viewOrder == "BU" ]]; then
        unWarpDir="-y"
      elif [[ $viewOrder == "TD" ]]; then
        unWarpDir="y"
      else
        unWarpDir=Unknown
      fi
    elif [[ $peDir == "RL" ]]; then
      if [[ $viewOrder == "BU" ]]; then
        unWarpDir="-x"
      elif [[ $viewOrder == "TD" ]]; then
        unWarpDir="x"
      else
        unWarpDir=Unknown
      fi
    else
      unWarpDir=Unknown
    fi

  else
    unWarpDir=Unknown
  fi

  ###############################


  ### sliceorder ################

    #SO=0=Sequential (seq)
    #SO=1=Interleaved (alt)


  #For now, multiband can only be interleaved
  if [[ $dcmSubSeries == "multiband_epi" || $dcmSubSeries == "multiband_mux_epi" ]]; then
    dcmSlicePattern="alternating"
    dcmTpattern="alt+z"
  else
    tmpsliceOrder=`gdcmdump -P $dcmPic | grep "SLICEORDER" | awk '{print $2}' | sed 's/\"//g'`
    if [[ $tmpsliceOrder -eq 0 ]]; then
      dcmSlicePattern="sequential"
      dcmTpattern="seq+z"
    else
      dcmSlicePattern="alternating"
      dcmTpattern="alt+z"
    fi
  fi

  ###############################
fi

#If the output base name is not specified, default to the Series Description in the DICOM header
if [[ $outBase == "" ]]; then
  outBase=${dcmSeries}
fi

#Spit out information about the DICOM series
if [[ $dcmSubSeries == "multiband_mux_epi" || $dcmSubSeries == "epiRT" || $dcmSubSeries == "multiband_epi" || $dcmSubSeries == "ia_stable_epiRTTopUp" ]]; then
  echo "SeriesNum=${dcmRelSeries}"
  echo ".Series=${dcmSeries}"
  echo "..TR=${dcmTR}"
  echo "...TE=${dcmTE}"
  echo "....${dcmNumVols} volume(s), ${dcmNumSlices} slices per volume"
  echo ".....SliceOrientation=${scanType}"
  echo ".....PhaseEncDir=${peDir}"
  echo ".....ViewOrder=${viewOrder}"
  echo ".....unWarpDir=${unWarpDir}"
  echo ".....dwellTime=${dwellTime}"
  echo ""
else
  echo "SeriesNum=${dcmRelSeries}"
  echo ".Series=${dcmSeries}"
  echo "..TR=${dcmTR}"
  echo "...TE=${dcmTE}"
  echo "....${dcmNumVols} volume(s), ${dcmNumSlices} slices per volume"
fi


#######################################################################
# Sort the list of files based on temporal position
#######################################################################

#Create a list of the DICOM files which will be probed for temporal position
#dcmList=`ls | grep \.dcm`
numDicom=${#dcmList[@]}

#Check for total number of DICOM files needed to convert
totalDicom=`expr $dcmNumSlices \* $dcmNumVols`
if [[ $totalDicom -ne $numDicom ]]; then
  echo "Error:  Expected $totalDicom images but found $numDicom"
  exit 1
fi

if [ -e $outDir/${outBase}_fileList.txt ]; then
  rm $outDir/${outBase}_fileList.txt
fi

if [ -e $outDir/${outBase}_tmpDCM.txt ]; then
  rm $outDir/${outBase}_tmpDCM.txt
fi

echo ""
echo "Sorting DICOM images"
echo ""

#make sure the file exists (just to be safe)
touch $outDir/${outBase}_tmpDCM.txt

#determine how much faster your computer can chug the data
num_processes=$(cpu_thread_num)

#0020,0013 : Slice Number
#0020,1041 : Location of slice in scanner space
#0021,105e : Acquisition time, per slice of an individual TR (this is not always reliable and often reports as "0")
  #printf '%s\n' "${dcmList[@]}" | xargs -n1 -P ${num_processes} -I {} bash -c 'dicom_hinfo -tag 0020,0013 0020,1041 0021,105e {}' > $outDir/${outBase}_tmpDCM.txt
printf '%s\n' "${dcmList[@]}" | xargs -n 1 -P ${num_processes} -I {} bash -c 'dicom_hinfo -tag 0020,0013 0020,1041 0021,105e 0018,1060 {}' >> $outDir/${outBase}_tmpDCM.txt

#sort by slice number
sort -k2 -n $outDir/${outBase}_tmpDCM.txt > $outDir/${outBase}_sortDCM.txt

#######################################################################
# Reorder the "bands" of data reconstructed for multiband
#######################################################################

#v1.0 multiband
#if [[ $dcmSubSeries == "multiband_mux_epi" ]]; then
  #i=25
  #while [[ $i -le $dcmTot ]];
    #do
    #j=`echo $i | awk '{print ($1+25)}'`
   # k=`echo $i | awk '{print ($1+50)}'`

    #Reorder the files and place in a temporary list
    #cat $outDir/${outBase}_sortDCM.txt | head -n+${k} | tail -n-25 >> $outDir/${outBase}_sortDCM2.txt
    #cat $outDir/${outBase}_sortDCM.txt | head -n+${j} | tail -n-25 >> $outDir/${outBase}_sortDCM2.txt
    #cat $outDir/${outBase}_sortDCM.txt | head -n+${i} | tail -n-25 >> $outDir/${outBase}_sortDCM2.txt

    #let i=i+75
  #done
  #mv $outDir/${outBase}_sortDCM2.txt $outDir/${outBase}_sortDCM.txt
#fi

#"multiband_epi" looks to be broken into four sections
      #1-15 (1-14) top
      #16-30 (15-28) middle-A
      #31-45 (29-42) middle-A
      #46-60 (43-56) bottom
#11/2017 multiband_mux_epi
  #New version has four bands, 14 slices per band, different stacking order
  #See notes well above (original "v1" mux is now commented out and not used)

#5/1/18 version.  52 slices per volume, down from 56 slices per volume
if [[ $dcmSubSeries == "multiband_mux_epi" ]]; then
  incr=`echo "$dcmNumSlices $dcmNumBands" | awk '{print $1/$2}'`
  i=$incr
  while [[ $i -le $dcmTot ]];
    do
    let j=$i+$incr
    let k=$j+$incr
    let l=$k+$incr

    #Reorder the files and place in a temporary list
    cat $outDir/${outBase}_sortDCM.txt | head -n+${i} | tail -n-${incr} >> $outDir/${outBase}_sortDCM2.txt
    cat $outDir/${outBase}_sortDCM.txt | head -n+${j} | tail -n-${incr} >> $outDir/${outBase}_sortDCM2.txt
    cat $outDir/${outBase}_sortDCM.txt | head -n+${k} | tail -n-${incr} >> $outDir/${outBase}_sortDCM2.txt
    cat $outDir/${outBase}_sortDCM.txt | head -n+${l} | tail -n-${incr} >> $outDir/${outBase}_sortDCM2.txt

    let i=i+$dcmNumSlices
  done
  mv $outDir/${outBase}_sortDCM2.txt $outDir/${outBase}_sortDCM.txt
fi

if [[ $dcmSubSeries == "multiband_epi" ]]; then
  i=15
  while [[ $i -le $dcmTot ]];
    do
    j=`echo $i | awk '{print ($1+15)}'`
    k=`echo $i | awk '{print ($1+30)}'`
    l=`echo $i | awk '{print ($1+45)}'`

    #Reorder the files and place in a temporary list
    cat $outDir/${outBase}_sortDCM.txt | head -n+${l} | tail -n-15 >> $outDir/${outBase}_sortDCM2.txt
    cat $outDir/${outBase}_sortDCM.txt | head -n+${k} | tail -n-15 >> $outDir/${outBase}_sortDCM2.txt
    cat $outDir/${outBase}_sortDCM.txt | head -n+${j} | tail -n-15 >> $outDir/${outBase}_sortDCM2.txt
    cat $outDir/${outBase}_sortDCM.txt | head -n+${i} | tail -n-15 >> $outDir/${outBase}_sortDCM2.txt

    let i=i+60
  done
  mv $outDir/${outBase}_sortDCM2.txt $outDir/${outBase}_sortDCM.txt
fi


#######################################################################
# Determine slice direction (ascending vs. descending)
#######################################################################

#Compare the first and last slice, in a given TR volume, using spatial location to determine:
  #Axial: IS/SI
  #Coronal: PA/AP
  #Sagittal: RL/LR

dcmSlicePos1=`cat $outDir/${outBase}_sortDCM.txt | head -n+1 | awk '{print $3}'`
dcmSlicePos2=`cat $outDir/${outBase}_sortDCM.txt | head -n+${dcmNumSlices} | tail -n-1 | awk '{print $3}'`

#Whereas, in a shell, 0=true and 1=false, we need to evaluate float values, therefore we have to use "bc."  But, bc logic follows that 0=false and 1=true
dcmSliceComp=`echo ${dcmSlicePos1}'<'${dcmSlicePos2} | bc -l`

if [[ $dcmSliceComp -eq 1 ]]; then
  dcmSliceDir="ascending"
else
  dcmSliceDir="descending"
fi


##########################################################################################################
#Determine order for AFNI's to3d
# alt+z = altplus    = alternating in the plus direction    (0  600  200  800  400)
# seq+z = seqplus    = sequential in the plus direction     (0  200  400  600  800)
  ##These are AFNI options, but GE won't allow
# alt+z2             = alternating, starting at slice #1    (400    0  600  200  800)
# alt-z = altminus   = alternating in the minus direction   (400  800  200  600    0)
# alt-z2             = alternating, starting at slice #nz-2 (800  200  600    0  400)
# seq-z = seqminus   = sequential in the minus direction    (800  600  400  200    0)
##########################################################################################################
  #This section only pertains to EPI images

#Echo out information about ordering of files
echo ""
echo "Slice data determined to have been collected ${dcmSliceDir}, ${dcmSlicePattern} order (${dcmTpattern})"
echo ""
echo "Converting data from DICOM to ${imageType}"

#Generate a list of just the DICOM files to convert
cat $outDir/${outBase}_sortDCM.txt | awk '{print $1}' > $outDir/${outBase}_fileList.txt

#For multiband, generate a list, for the first band, of the slice-timing needed for adjustment, adjust for the remaining bnads
if [[ $dcmSubSeries == "multiband_mux_epi" || $dcmSubSeries == "multiband_epi" ]]; then
  bandSlices=`echo "$dcmNumSlices $dcmNumBands" | awk '{print $1/$2}'`
  cat $outDir/${outBase}_sortDCM.txt | head -n+${bandSlices} | awk '{print $5}' > $outDir/${outBase}_tmptimeList.txt

  #Need to adjust so that slice1 has time 0, not what is reported by ACQ Trigger Time (e.g. 76 ms should be 0), adjust all other slices by this time (slice2, 615 ms becomes 539)
  offsetTime=`cat $outDir/${outBase}_tmptimeList.txt | head -n+1 | tail -n-1`

  i=1
  while [[ $i -le $bandSlices ]];
    do
    if [[ $i -eq 1 ]]; then
      echo "0" >> $outDir/${outBase}_tmptimeList2.txt
    else
      cat $outDir/${outBase}_tmptimeList.txt | head -n+${i} | tail -n-1 | awk -v var=${offsetTime} '{print ($1-var)/1000}' >> $outDir/${outBase}_tmptimeList2.txt
    fi
  let i=i+1
  done

  #Now adjust for the multiple bands (already created one band, need to repeat for remainder)
  leftoverBands=`echo $dcmNumBands | awk '{print $1-1}'`
  cp $outDir/${outBase}_tmptimeList2.txt $outDir/${outBase}_timeList.txt
  i=1
  while [[ $i -le $leftoverBands ]];
    do
    cat $outDir/${outBase}_tmptimeList2.txt >> $outDir/${outBase}_timeList.txt
  let i=i+1
  done

  #Remove tmp timelist files
  rm $outDir/${outBase}_tmptimeList.txt $outDir/${outBase}_tmptimeList2.txt
fi

#Run the conversion
if [[ $dcmSubSeries == "multiband_mux_epi" || $dcmSubSeries == "epiRT" || $dcmSubSeries == "multiband_epi" || $dcmSubSeries == "ia_stable_epiRTTopUp" ]]; then
  if [[ $dcmNumVols -eq 1 ]]; then
    dcmPicOne=`cat $outDir/${outBase}_fileList.txt | head -n+1`
    mri_convert --in_type dicom --out_type nii $dcmPicOne $outDir/${outBase}.nii.gz
  else

    if [[ $dcmVerbose == 0 ]]; then
      comto3d="to3d -session $outDir -prefix ${outBase}${imageExt} -time:zt $dcmNumSlices $dcmNumVols $dcmTR $dcmTpattern -@ < $outDir/${outBase}_fileList.txt"
      to3d -session $outDir -prefix ${outBase}${imageExt} -time:zt $dcmNumSlices $dcmNumVols $dcmTR $dcmTpattern -@ < $outDir/${outBase}_fileList.txt
    else
      comto3d="to3d -session $outDir -save_outliers $outDir/${outBase}_outliers.1d -prefix ${outBase}${imageExt} -time:zt $dcmNumSlices $dcmNumVols $dcmTR $dcmTpattern -@ < $outDir/${outBase}_fileList.txt"
      to3d -session $outDir -save_outliers $outDir/${outBase}_outliers.1d -prefix ${outBase}${imageExt} -time:zt $dcmNumSlices $dcmNumVols $dcmTR $dcmTpattern -@ < $outDir/${outBase}_fileList.txt
      outBase2=`echo $outBase | sed -e 's/_/ /g'`
      1dplot -png $outDir/${outBase}_outliers.png -title "'${outBase2}' Outliers" -xlabel "TR" -ylabel "Intensity" -one $outDir/${outBase}_outliers.1d
    fi
    #For multiband, have to adjust slice-timing as "bands" aren't accounted for with to3d
    if [[ $dcmSubSeries == "multiband_mux_epi" || $dcmSubSeries == "multiband_epi" ]]; then
      com3dTshift="3dTshift -TR $dcmTR -prefix ${outDir}/${outBase}_tmp${imageExt} -quintic -tpattern @${outDir}/${outBase}_timeList.txt -ignore 3 ${outDir}/${outBase}${imageExt}"
      3dTshift -TR $dcmTR -prefix ${outDir}/${outBase}_tmp${imageExt} -quintic -tpattern @${outDir}/${outBase}_timeList.txt -ignore 3 ${outDir}/${outBase}${imageExt}
      mv ${outBase}_tmp${imageExt} ${outBase}${imageExt}
      #For now, saving timing file
        #rm $outDir/${outBase}_timeList.txt
    fi
  fi
#Just a localizer, new sequence is wonky and AFNI gripes, just brute force it into a 3D volume via mri_convert
#elif [[ $dcmSubSeries == "ssfse" ]]; then
#  comto3d="to3d -session $outDir -prefix ${outBase}${imageExt} -time:zt $dcmNumSlices $dcmNumVols $dcmTR $dcmTpattern -@ < $outDir/${outBase}_fileList.txt"
#  to3d -session $outDir -prefix ${outBase}${imageExt} -time:zt $dcmNumSlices $dcmNumVols $dcmTR $dcmTpattern -@ < $outDir/${outBase}_fileList.txt
else
  #Source the first DICOM in the series
  dcmPicOne=`cat $outDir/${outBase}_fileList.txt | head -n+1`
  #convert via mri_convert
  mri_convert --in_type dicom --out_type nii $dcmPicOne $outDir/${outBase}.nii.gz

  #Issue where fieldMap is RAS but comes out as RAI
    #https://afni.nimh.nih.gov/afni/community/board/read.php?1,79549,79553#msg-79553
  if [[ $dcmSubSeries == "ia_stable_B0map" ]]; then
    3drefit -orient RAS $outDir/${outBase}.nii.gz
  fi

# Keeping afni info in nifti header caused headaches
# see: https://afni.nimh.nih.gov/afni/community/board/read.php?1,146771,146771#msg-146771
# removing afni info from nifti files
#check for afni extensions
ext_check=`nifti_tool -disp_ext -infiles $outDir/${outBase}.nii.gz -debug 1 | awk -F"= " '{print $2}'`
echo "extension check: ${ext_check}"
if [[ "${ext_check}" > "0" ]]; then
  mv $outDir/${outBase}.nii.gz $outDir/tmp_${outBase}.nii.gz
  nifti_tool -rm_ext ALL -infiles $outDir/tmp_${outBase}.nii.gz -prefix $outDir/${outBase}.nii.gz
  rm $outDir/tmp_${outBase}.nii.gz
fi

  #If AFNI HEAD/BRIK is desired, convert and remove NIfTI file
  if [[ $imageType == "AFNI" ]]; then
    3dcopy $outDir/${outBase}.nii.gz $outDir/${outBase}
    rm $outDir/${outBase}.nii.gz
  fi
fi

#Check the native orientation, for logging
native_orient $outDir/${outBase}${imageExt}

#Reorient file to RPI (only if selected via the "-r" flag)
  #This will be skipped if the file is a localizer
if [[ $reorientFlag == 1 ]]; then
  if [[ $dcmSubSeries != "ssfse" ]]; then
    if [[ $imageType == "NIfTI" ]]; then
      RPI_orient $outDir/${outBase}${imageExt}
    else
      3dresample -orient rpi -rmode Cu -prefix $outDir/tmp -inset  $outDir/${outBase}+orig
      mv $outDir/tmp+orig.HEAD $outDir/${outBase}+orig.HEAD
      mv $outDir/tmp+orig.BRIK $outDir/${outBase}+orig.BRIK
    fi
  fi
fi

#If multiband EPI, strip off the first volume (average calibration scan), or adjust for new sequences (e.g. mux)
  #Based on voxel intensities for mux sequence, first 3 TRs look to be adequate for removal
if [[ $dcmSubSeries == "multiband_mux_epi" || $dcmSubSeries == "multiband_epi" ]]; then
  if [[ $imageType == "NIfTI" ]]; then
    3dTcat -session $outDir -prefix tmp.nii.gz $outDir/${outBase}${imageExt}'[4..$]'
    3dTcat -session $outDir -prefix tmp1.nii.gz $outDir/${outBase}${imageExt}'[0..3]'
    mv $outDir/tmp.nii.gz $outDir/${outBase}${imageExt}
    mv $outDir/tmp1.nii.gz $outDir/disDacq${imageExt}
  else
    3dTcat -session $outDir -prefix tmp $outDir/${outBase}+orig'[4..$]'
    3dTcat -session $outDir -prefix tmp1 $outDir/${outBase}+orig'[0..3]'
    mv $outDir/tmp+orig.HEAD $outDir/${outBase}+orig.HEAD
    mv $outDir/tmp+orig.BRIK $outDir/${outBase}+orig.BRIK
    mv $outDir/tmp1+orig.HEAD $outDir/disDacq+orig.HEAD
    mv $outDir/tmp1+orig.BRIK $outDir/disDacq+orig.BRIK
  fi
fi

#####################################################################################
# Special Processing for fieldmaps to make them BIDS compatible (GE)
#####################################################################################
if [[ "${dcmSubSeries}" == "ia_stable_B0map" ]]; then
  mv $outDir/${outBase}_RPI${imageExt} $outDir/${outBase}_master${imageExt}
  fslroi $outDir/${outBase}_master${imageExt} $outDir/${outBase}${imageExt} 0 1
  fslroi $outDir/${outBase}_master${imageExt} $outDir/${outBase}_magnitude${imageExt} 1 1
  rm $outDir/${outBase}_master${imageExt}
  # fieldmap should be a float32 not an int16
  fslmaths $outDir/${outBase}${imageExt} $outDir/${outBase}${imageExt} -odt float
fi

#####################################################################################
# Create a log to save information about conversion
#####################################################################################

echo ""
echo "Logging DICOM information"
echo ""

#Retain information about the DICOM data
echo "${dcmSeries}/${dcmSubSeries}" > $outDir/${outBase}_info.txt
echo "" >> $outDir/${outBase}_info.txt
echo "scannerType=${scanMan}/${scanModel}" >> $outDir/${outBase}_info.txt
echo "magnetStrength=${magSize}T" >> $outDir/${outBase}_info.txt
echo "SeriesNum=${dcmRelSeries}" >> $outDir/${outBase}_info.txt
echo "TR=${dcmTR}" >> $outDir/${outBase}_info.txt
echo "TE=${dcmTE}" >> $outDir/${outBase}_info.txt
echo "sliceThickness=${sliceThickness}" >> $outDir/${outBase}_info.txt
echo "flipAngle=${flipAng}" >> $outDir/${outBase}_info.txt
echo "echoSpacing=${echoSpacing}" >> $outDir/${outBase}_info.txt
echo "matrixSize=${matrixSize}" >> $outDir/${outBase}_info.txt
echo "${dcmNumVols} volumes, ${dcmNumSlices} slices per volume" >> $outDir/${outBase}_info.txt

#Print out relevent information for EPI data
if [[ $dcmSubSeries == "multiband_mux_epi" || $dcmSubSeries == "epiRT" || $dcmSubSeries == "multiband_epi" || $dcmSubSeries == "ia_stable_epiRTTopUp" ]]; then
  echo "SliceOrientation=${scanType}" >> $outDir/${outBase}_info.txt
  echo "PhaseEncodingDirection=${peDir}" >> $outDir/${outBase}_info.txt
  echo "ViewOrder=${viewOrder}" >> $outDir/${outBase}_info.txt
  echo "unWarpDir=${unWarpDir}" >> $outDir/${outBase}_info.txt
  echo "dwellTime=${dwellTime}" >> $outDir/${outBase}_info.txt
  echo "SliceOrder=${dcmSliceDir}" >> $outDir/${outBase}_info.txt
  echo "dcmSlicePattern=${dcmSlicePattern}" >> $outDir/${outBase}_info.txt
  echo "dcmTpattern=${dcmTpattern}" >> $outDir/${outBase}_info.txt
  echo "NativeOrientation=${native_orient}" >> $outDir/${outBase}_info.txt
  echo "" >> $outDir/${outBase}_info.txt
  echo "to3d command:" >> $outDir/${outBase}_info.txt
  echo "  ${comto3d}" >> $outDir/${outBase}_info.txt
  if [[ $dcmSubSeries == "multiband_mux_epi" || $dcmSubSeries == "multiband_epi" ]]; then
    echo "  ${com3dTshift}" >> $outDir/${outBase}_info.txt
  fi
  echo "" >> $outDir/${outBase}_info.txt
  echo "" >> $outDir/${outBase}_info.txt
  echo "" >> $outDir/${outBase}_info.txt
  echo "DICOM Image : Slice Scanner_Postion Slice_Acquisition_Time Slice_Trigger_Time" >> $outDir/${outBase}_info.txt
  echo "-------------------------------------------------------------------------------" >> $outDir/${outBase}_info.txt
  cat $outDir/${outBase}_sortDCM.txt >> $outDir/${outBase}_info.txt
fi


if [ ! -z ${BIDS} ]; then
  # making a JSON file
  # Assuming TR and TE are in milliseconds
  TR=$(echo "scale=3; ${dcmTR}/1000" | bc | awk '{printf "%f", $0}')
  TE=$(echo "scale=10; ${dcmTE}/1000" | bc | awk '{printf "%f", $0}')

  # change phase encoding direction to i,j,k, and use rev if warp is negative
  BIDS_unWarpDir=$(echo ${unWarpDir} | sed -e "s|x|i|" -e "s|y|j|" -e "s|z|k|" | rev)
  # total readout time = 1/pixelbandwidth
  # tag 0018,0095
  PixelBandWidth=$(dicom_hinfo -no_name -tag "0018,0095" $dcmPic)
  # Email Joel about total Readout time
  # TotalReadOutTime=$(echo "1/${PixelBandWidth}" | bc -l)
  # The EffectiveEchoSpacing refers to dwell time for processing
   if [[ $dcmSubSeries == "multiband_mux_epi" || $dcmSubSeries == "epiRT" || $dcmSubSeries == "multiband_epi" || $dcmSubSeries == "ia_stable_epiRTTopUp" ]]; then
      echo "{
      \"Manufacturer\": \"${scanMan}\",
      \"ManufacturersModelName\": \"${scanModel}\",
      \"RawImage\": false,
      \"SeriesDescription\": \"${dcmSeries}\",
      \"MagneticFieldStrength\": ${magSize},
      \"FlipAngle\": ${flipAng},
      \"EchoTime\": ${TE},
      \"EffectiveEchoSpacing\": ${dwellTime},
      \"PhaseEncodingDirection\": \"${BIDS_unWarpDir}\",
      \"RepetitionTime\": ${TR},
      \"ConversionSoftware\": \"GE_dcm_to_nii.sh\"
}" > $outDir/${outBase}.json
  elif [[ "${dcmSubSeries}" == "ia_stable_B0map" ]]; then
    echo "{
     \"Manufacturer\": \"${scanMan}\",
     \"ManufacturersModelName\": \"${scanModel}\",
     \"RawImage\": false,
     \"SeriesDescription\": \"${dcmSeries}\",
     \"MagneticFieldStrength\": ${magSize},
     \"Units\": \"Hz\",
     \"ConversionSoftware\": \"GE_dcm_to_nii.sh\"
}" > ${outDir}/${outBase}.json
   else
     echo "{
      \"Manufacturer\": \"${scanMan}\",
      \"ManufacturersModelName\": \"${scanModel}\",
      \"RawImage\": false,
      \"SeriesDescription\": \"${dcmSeries}\",
      \"MagneticFieldStrength\": ${magSize},
      \"FlipAngle\": ${flipAng},
      \"EchoTime\": ${TE},
      \"RepetitionTime\": ${TR},
      \"ConversionSoftware\": \"GE_dcm_to_nii.sh\"
}" > ${outDir}/${outBase}.json
  fi
fi

#Cleanup
rm $outDir/${outBase}_sortDCM.txt
rm $outDir/${outBase}_tmpDCM*.txt

popd > /dev/null
popd > /dev/null

