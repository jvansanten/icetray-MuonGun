#!/bin/sh

# -------------------------------------------------------
# GridFTP wrapper script, based on an example by Heath
# Skarlupa: https://wiki.icecube.wisc.edu/index.php/GZK_Cluster
# -------------------------------------------------------

uid=51509

function print_usage () {
  echo "  options:"
  echo "    -t TARBALL  IceTray tarball to unpack"
  echo
  echo "    -d DEST     Destination directory: files will be copied to "
  echo "                gsiftp://data.icecube.wisc.edu/DEST/BASENAME"
  echo
  echo "    -i          Switch to input mode. All following arguments will be interpreted"
  echo "                as files to copy from gsiftp://data.icecube.wisc.edu"
  echo
  echo "    -o          Switch to output mode. All following arguments will be interpreted"
  echo "                as files to copy back to gsiftp://data.icecube.wisc.edu"
  echo
  echo "    --          End argument parsing. All following arguments will be"
  echo "                run as a script."
  exit 1
}

if [[ "$#" == "0" ]]; then
  print_usage
fi

next="input"
while (( "$#" )); do
  # break out if we've seen "--"
  if [[ $next == "args" ]]; then
    break
  fi
  case $1 in
    "-d")
      shift;
      dest=$1
      shift;;
    "-t")
      shift;
      tarball=$1
      shift;;
    "-i")
      next="input";
      shift;;
    "-o")
      next="output";
      shift;;
    "--")
      next="args";
      shift;;
    *)
      case $next in 
        "input")
          inputs="$inputs $1";
          shift;;
        "output")
          outputs="$outputs $1";
          shift;;
        "args")
          break;;
      esac
  esac
done

if [[ -z $dest ]]; then
  echo "You must specify a destination!"
  print_usage
fi

if [[ -z $tarball ]]; then
  echo "You must specify a tarball!"
  print_usage
fi

# Try to figure out where the heck JAVA_HOME is
for candidate in /home/icecube/tools/j2sdk1.4.2-x86_64 /data2/icecube/tools/j2sdk1.4.2-x86_64 /usr/java/j2sdk1.4.2_12 /home/icecube/i3tools/j2sdk1.4.2_12; do
  if [[ -d $candidate ]]; then
    export JAVA_HOME=$candidate
    break
  fi
done

if [[ -z $JAVA_HOME ]]; then
  echo "Can't find JAVA_HOME!"
  exit 1
fi

# Unpack certificates
tar xvzf doegrids.tar.gz
# Location of CA Certs Directory
X509_CERT_DIR=${_CONDOR_SCRATCH_DIR}/doegrids
export X509_CERT_DIR
# Location of my user proxy
X509_USER_PROXY=${_CONDOR_SCRATCH_DIR}/x509up_u${uid}
export X509_USER_PROXY

# Transfer inputs
for infile in $inputs $tarball; do
  globus-url-copy -p 8 -vb gsiftp://data.icecube.wisc.edu/$infile file://${_CONDOR_SCRATCH_DIR}/`basename $infile`
done

# Unpack the tarball
tar xvzf `basename $tarball`

export I3_BUILD=${_CONDOR_SCRATCH_DIR}/`basename $tarball .tar.gz`

# Run the script in an env-shell
$I3_BUILD/env-shell.sh $@

# Copy outputs back
for outfile in $outputs; do
  globus-url-copy -p 8 -vb file://${_CONDOR_SCRATCH_DIR}/$outfile gsiftp://data.icecube.wisc.edu/$dest/$outfile
done
