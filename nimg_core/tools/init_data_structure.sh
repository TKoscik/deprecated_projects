#! /bin/bash

# Help ------------------------------------------------------------------------
function Usage {
    cat <<USAGE

`basename $0` initializes a data structure for a neuroimaging project.

Usage:
`basename $0`  -r researcher_directory
                        -p project_name
                        -g group_permission
                        -h <help>

Example:
  bash $0 -r /Shared/koscikt_scratch -p example_data

Arguments:
  -r researcher_directory  The full directory listing in which to create the
                           project data structure.
  -p project_name          A unique name for the imaging project, must not
                           exist.
  -g group_permission      A string indicating the group to assign directory
                           ownership and access
  -h help

USAGE
    exit 1
}

# Parse input operators -------------------------------------------------------
while getopts "r:p:g:h" option
do
  case "${option}"
  in
    r) # researcher_directory
      researcherRoot=${OPTARG}
      ;;
    p) # project_name
      projectName=${OPTARG}
      ;;
    g) # group permissions
      groupOwn=${OPTARG}
      ;;
    h) # help
      Usage >&2
      exit 0
      ;;
    *) # unknown options
      echo "ERROR: Unrecognized option -$OPT $OPTARG"
      exit 1
      ;;
  esac
done

# Create Directory tree -------------------------------------------------------
if [ -d "${researcherRoot}" ]; then
  if [ ! -d "${researcherRoot}/${projectName}" ]; then
    mkdir ${researcherRoot}/${projectName}
    mkdir ${researcherRoot}/${projectName}/dicom
    mkdir ${researcherRoot}/${projectName}/nifti
    mkdir ${researcherRoot}/${projectName}/deriv
    mkdir ${researcherRoot}/${projectName}/scripts
    mkdir ${researcherRoot}/${projectName}/summary
    chgrp -R ${groupOwn} ${researcherRoot}/${projectName}
    chmod -R g+rw ${researcherRoot}/${projectName}
  else
    echo "ERROR: ${researcherRoot}/${projectName} already exists"
  fi
else
  echo "ERROR: ${researcherRoot} does not exist"
fi
