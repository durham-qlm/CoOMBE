# Examples

This folder contains the files listed below. Using these programs normally requires the modules contained in the ldbl, obe and mbe files, as well as (except for the bespoke programs) the driveall program. The general_settings module is also necessary; a copy of this module is provided with each example, with the number of states adapted to the requirements of this example. Besides the files listed below, each folder also contains the **Makefile** file, **Dockerfile** file and **.sh** shell script used when running the program with Podman as explained in Appendix G of the article. The Makefile file and the shell scipt can also be used for compiling and running the codes outside Podman.  

## Documentation

- **examples.pdf** : Further examples of the use of the driveall program for running the code.

## Files related to the example given in Appendix E of the article (folder ex_example)

- **general_settings_3st.f90**  : The general_settings module for a 3-state model.
- **example.f90**       : The bespoke example program (the same example is also given in Section 3.7 of the User Manual).
- **example_c.dat**     : The controlparams file for this example.
- **example_k.dat**     : The keyparams file for this example.


## Files related to the example given in Appendix F of the article (folder ex_propagation)

- **general_settings_3st.f90** : The general_settings module for a 3-state model.
- **prop.f90** : A bespoke program, not based on the driveall program, doing the same propagation calculation as that described in this appendix.
- **prop_c.dat** : The controlparams file for this example.
- **prop_k.dat** : The keyparams file for this example.
- **appliedfields.dat** : The tdamps_in file for this example.
- **outamplitudes.dat** : File produced by the driveall program in this example.
- **fig_prop.py** : The python program used for plotting the figure.
- **figure_prop.pdf** : The output of fig_prop.py.



## Files related to the first additional example, section 2 of examples.pdf (folder ex_4st_ladder)

- **general_settings_4st.f90**  : The general_settings module for a 4-state model.
- **ladder_c.dat**      : The controlparams file for this example. 
- **ladder_k.dat**      : The keyparams file for this example.
- **ladder_rho21.dat**  : File produced by the driveall program in this example.
- **fig_ladder.py**     : The python program used for plotting Fig. 1 of example.pdf.
- **figure_ladder.pdf** : The output of fig_ladder.py.



## Files related to the second additional example, section 3 of examples.pdf (folder ex_abs_D1_Rb)

- **general_settings_24st.f90** : The general_settings module for a 24-state model.
- **Rb85D1_c.dat**      : The controlparams file for this example.
- **Rb85D1_d.dat**      : The defaultdata file for this example.
- **Rb85D1_k.dat**      : The keyparams file for this example.
- **chi_noweakprb**     : File produced by the driveall program in this example.
- **chi_weakprb_2**     : File produced by the driveall program in this example.
- **chi_weakprb_3**     : File produced by the driveall program in this example.
- **fig_Rb85D1.py**     : The python program used for plotting Fig. 2 of example.pdf.
- **figure_Rb85D1.pdf** : The output of fig_Rb85D1.py.



## Files related to the third additional example, section 4 of examples.pdf (folder ex_rabi_osc)

- **general_settings_2st.f90**  : The general_settings module for a 2-state model.
- **timedep_c.dat**      : The controlparams file for this example.
- **timedep_k.dat**      : The keyparams file for this example.
- **timedep_rho11.dat**  : File produced by the driveall program in this example.
- **timedep_rho12.dat**  : File produced by the driveall program in this example.
- **fig_timedep.py**     : The python program used for plotting Fig. 3 of example.pdf.
- **figure_timedep.pdf** : The output of fig_timedep.py.
