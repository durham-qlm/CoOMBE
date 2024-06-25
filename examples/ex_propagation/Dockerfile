#create base podman image using arch linux
FROM archlinux/archlinux:base-devel
RUN pacman -Syyu --noconfirm
RUN pacman-db-upgrade
#install fortran compiler and libraries
RUN pacman -S --noconfirm make sudo gcc-fortran blas lapack
#create an output directory
RUN mkdir -p output
#copy program source code, and input files to podman image
COPY Makefile *.f90 *.sh *.dat /CoOMBE/
#change working directory to podman image
WORKDIR /CoOMBE/
#compile program using Makefile
RUN make
#Run program, then copy output file from podman image to local output directory
CMD ./run_propagation.sh 
