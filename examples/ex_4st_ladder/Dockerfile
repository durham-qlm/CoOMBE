#create base podman image using arch linux
FROM archlinux/archlinux:base-devel
RUN pacman -Syyu --noconfirm
RUN pacman-db-upgrade
#install fortran compiler and libraries
RUN pacman -S --noconfirm make sudo gcc-fortran blas lapack
#create home directory, copy program source code and input files to image
RUN mkdir -p /home/coombe
COPY Makefile *.f90 /home/coombe/ 
COPY *.sh *.dat /home/coombe/
WORKDIR /home/coombe/
#Run program, then copy output file from podman image to local output directory
CMD ./run_4st_ladder.sh 
