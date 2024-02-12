#!/bin/bash

#SCRIPT_PATH="${BASH_SOURCE:-$0}"
#ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
##echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
#ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
##echo "Value of ABS_DIRECTORY: ${ABS_DIRECTORY}"

echo "Installing pre-requisite tools.."
bash basic_setup_ubuntu.sh
echo "Done"

echo "Downloading and setting up miniconda..."
if [ ! -e "Miniconda3-py39_23.3.1-0-Linux-x86_64.sh" ]
then
    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.3.1-0-Linux-x86_64.sh
fi

bash ./Miniconda3-py39_23.3.1-0-Linux-x86_64.sh -b -p ./miniconda3
echo "Downloading and setting up miniconda...DONE"

echo "Seeting up conda env named with given argument"
miniconda3/bin/conda env create --name distbwa -f environment.yml
echo "Seeting up conda env named new_env...DONE"

echo "Activating conda env..."
source miniconda3/bin/activate distbwa
echo "localhost" > hostfile

## build tools
WDIR=../../
EXEDIR=`pwd`

# compile htslib
# cd ${WDIR}/applications/htslib
# autoreconf -i  # Build the configure script and install files it uses
# ./configure    # Optional but recommended, for choosing extra functionality
# make
#make install   #uncomment this for installation

# compile samtools
# cd ${WDIR}/applications/samtools
# autoheader
# autoconf -Wno-syntax
# chmod 775 configure
# ./configure           # Needed for choosing optional functionality
# make
# saminstall="SUCESS"
# if [ -e "${WDIR}/applications/samtools/samtools" ]; then
#     echo "SAMTools build successful"
# else
#     saminstall="FAILED"
#     echo "Error!! SAMTools build failed"
# fi

# if [ "$?" == "0" ]
# then
#     echo "Samtools installed successfully"
# else
#     echo "Samtools installation failed"
# fi
#make install         #uncomment this for installation

# compile starsolo
echo "Build STAR" 
cd ${WDIR}/applications
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11a.tar.gz
tar -xvf 2.7.11a.tar.gz
rm 2.7.11a.tar.gz
mv STAR-2.7.11a/ STAR
cd STAR/bin/Linux_x86_64_static
chmod +x ./STAR

# compile samtools
echo "Build SAMTOOLS"
echo "EXEDIR"
cd $EXEDIR
pwd
cd ${WDIR}/applications
pwd
wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
tar -xf samtools-1.11.tar.bz2
cd samtools-1.11
./configure
make
make install
cd ..
mv samtools-1.11/ samtools
saminstall="SUCESS"
if [ -e "${WDIR}/applications/samtools/samtools" ]; then
    echo "SAMTools build successful"
else
    saminstall="FAILED"
    echo "Error!! SAMTools build failed"
fi

# compile htslib
echo "EXEDIR"
cd $EXEDIR
cd ${WDIR}/applications/htslib
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make

# fastqprocess
cd $EXEDIR
git clone --recursive https://github.com/broadinstitute/warp-tools.git -b develop
cd warp-tools/tools/fastqpreprocessing/
./fetch_and_make_dep_libs.sh && make
## make -j

if [ "$?" == "0" ]
then
    echo "fqprocess installed successfully"
else
    echo "fqprocess installation failed"
fi

echo "star compilation is "$starinstall
echo "samtools compilation is "$saminstall

echo "Compelete installation done."
