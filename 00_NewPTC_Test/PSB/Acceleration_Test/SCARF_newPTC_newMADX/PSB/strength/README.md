# PS Booster Optics Repository

## strength directory

Contains the magnetic element strength definition files for the MAD-X program

List of the official strength files, to be used with the corresponding scripts in the madx directory

* psb_orbit.str : the main file downloaded from the LayoutDB
* psb_injection.str : the modifications to the psb_orbit.str file relevant for the injection. 
* psb_extraction.str : the modificaitons to the psb_orbit.str file relevant for the extraction.

**Note**: all files should be loaded AFTER the psb_orbit.str file in the madx scripts.
