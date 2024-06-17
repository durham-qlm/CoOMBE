FROM archlinux/archlinux:base-devel

RUN pacman -Syyu --noconfirm
RUN pacman-db-upgrade
RUN pacman -S --noconfirm make sudo gcc-fortran blas lapack
RUN mkdir -p output

COPY Makefile *.f90 *.sh *.dat /honours/
WORKDIR /honours/

RUN make

CMD ./run_4st.sh ; cp /honours/ladder_rho21.dat /output/ladder_rho21.dat
#CMD ["ladder_rho21.dat"]
