# ANN_Classification_of_M12_resonant_argument_images
A repository for the use of artificial neural networks to automatically classify images of the resonant argument of the M1:2 mean-motion resonance.
The codes were developed for a Linux Ubuntu operating system, with gfortran, python3.8.5 and tensorflow installed.  They may be portability issues for other platforms.  In the main branch there is a UNIX script (script_analysis) that will automatically perform the analysis of images of resonant arguments for simulations in the RUN_00/01 directory. 

A sample of 99 images of resonant arguments and their labels, used to train the neural networks, are available in the branch ALL_PNG_varpic.  The full database of 5700 images and labels used in Carruba et al. (2021), ..., is available at ... .  Several codes used for the numerical simulations and their analysis are available
in the branch CODES.  swift_bs.f is a Burlisch-Stoer integrator from the SWIFT package (Levison and Duncan 1994, Icarus 108, 18, available at https://www.boulder.swri.edu/~hal/swift.html).  follow_all_filtered_plan_el.f will transform the outcome of the numerical simulation from a binary file to a ascii format. res_arg_m12.f will compute the resonant arguments for the M1:2 resonance for 50 test particles in the simulation (it needs a file saturn to be in the same 
directory). plot_rez_id.py will plot the time behavior of the resonant arguments.  image_class.py	will classify the images based on training obtained from images and labels at the image_class.py branch.  Finally, an example of a simulation of 50 asteroids with all the necessary input and output files is available in the RUN_00/01 branch.

In order to run the model, users will only have to run the script, after changing its ṕermissions, using the following commands:

chmod+x script_analysis

./script_analysis
