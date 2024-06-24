This folder contains the files listed below. Using these programs requires the
ldbl, obe and (in most cases) mbe modules, as well as (except for the bespoke
programs) the driveall program.

1. Documentation
----------------

examples.pdf              : Further examples of the use of the driveall program
                            for running the code.



2. Files related to the example given in Appendix E of the article
------------------------------------------------------------------

compile_3st       : Example of compilation statement for this case.
general_settings_3st.f90  : The general_settings module for a 3-state model.
example.f90       : The bespoke example program (the same example is also
                    given in Section 3.7 of the User Manual).
example_c.dat     : The controlparams file for this example.
example_k.dat     : The keyparams file for this example.



3. Files related to the example given in Appendix F of the article
------------------------------------------------------------------

compile_3st       : Example of compilation statement for this case.
general_settings_3st.f90  : The general_settings module for a 3-state model.
prop.f90          : A bespoke program, not based on the driveall program,
                    doing the same propagation calculation as that described
                    in this appendix.
prop_c.dat        : The controlparams file for this example.
prop_k.dat        : The keyparams file for this example.
appliedfields.dat : The tdamps_in file for this example.
outamplitudes.dat : File produced by the driveall program in this example.
fig_prop.py       : The python program used for plotting the figure.
figure_prop.pdf   : The output of fig_prop.py.



4. Files related to the first additional example (Section 2 of examples.pdf)
----------------------------------------------------------------------------

compile_4st       : Example of compilation statement for this case.
general_settings_4st.f90  : The general_settings module for a 4-state model.
ladder_c.dat      : The controlparams file for this example. 
ladder_k.dat      : The keyparams file for this example.
ladder_rho21.dat  : File produced by the driveall program in this example.
fig_ladder.py     : The python program used for plotting Fig. 1 of example.pdf.
figure_ladder.pdf : The output of fig_ladder.py.



5. Files related to the second additional example (Section 3 of examples.pdf)
-----------------------------------------------------------------------------

compile_24st      : Example of compilation statement for this case.
general_settings_24st.f90 : The general_settings module for a 24-state model.
Rb85D1_c.dat      : The controlparams file for this example.
Rb85D1_d.dat      : The defaultdata file for this example.
Rb85D1_k.dat      : The keyparams file for this example.
chi_noweakprb     : File produced by the driveall program in this example.
chi_weakprb_2     : File produced by the driveall program in this example.
chi_weakprb_3     : File produced by the driveall program in this example.
fig_Rb85D1.py     : The python program used for plotting Fig. 2 of example.pdf.
figure_Rb85D1.pdf : The output of fig_Rb85D1.py.



6. Files related to the third additional example (Section 4 of examples.pdf)
----------------------------------------------------------------------------

compile_2st        : Example of compilation statement for this case.
general_settings_2st.f90  : The general_settings module for a 2-state model.
timedep_c.dat      : The controlparams file for this example.
timedep_k.dat      : The keyparams file for this example.
timedep_rho11.dat  : File produced by the driveall program in this example.
timedep_rho12.dat  : File produced by the driveall program in this example.
fig_timedep.py     : The python program used for plotting Fig. 3 of example.pdf.
figure_timedep.pdf : The output of fig_timedep.py.
