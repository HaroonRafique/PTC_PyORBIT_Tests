# PS Booster Optics Repository

## madx directory

Contains the command scripts to run the MAD-X program


The extension for the MADX files is \<\>.madx and not \<\>.mad.

The reason is that when the files are read via the web, the .mad extension is taken as Microsoft Access module and not as
a text file. By using .madx extension, the files are interpreted as text files on the web.

Excecute the madx files from the "madx" directory, e.g.:
> >madx  < psb_extraction.madx

### MAD-X script files
List of official madx scripts:

* **psb_injection.madx**                : injection optics
* **psb_extraction_QH4.17_QV4.23.madx** : extraction optics with {4.17,4.23} tune
* **psb_orbit.madx**                    : orbit evaluations
* **psb_survey.madx**                   : to generate the survey files  

### Additional Information on MAD-X files

For the PSB there are two injection bumps: The chicane and the movable bump. The movable bump starts with a large value and is then decreased to zero.

The PS Booster ring consists of 16 identical periods, apart from the equipment in the straight sections. The BOOSTER has a radius of 25m and thus a circonference of 157.08m

The injected beam from the linac have Ekin=50 MeV proton beam with p~300MeV/c.
 (or to be exacrt 311 MeV/c according to the [PAC05 report](http://accelconf.web.cern.ch/accelconf/p05/PAPERS/TPAT054.PDF))
corresponding to E=0.988471 GeV for E0=0.938272 GeV. (Reminder : the rest mass of a proton is 0.938272 GeV/c2)

 The kinetic energy of the extracted beam is Ekin=1.4 GeV and p=2.141766 GeV/c, estimated from the formula
 > p=SQRT(Ekin^2+2\*E0\*Ekin)/c

 The beam is extracted to the PS. In the PS, it is injected in SS42.

 * Revolution time at extraction:  T,revolution, ext = 572.79 ns (nano  seconds)
 * Revolution time at injection:   T,revolution, inj =   1.6  us (micro seconds)

The BOOSTER today has 2 bunches in each ring, i.e. 8 bunches in the PS.

 * Injection time: C=275 ms, measure time C=290 ms

Working points:
 * Injection - high intensity          : QX=4.28,  QY=4.55
 * Injection - low  intensity e.g. LHC : QX=4.28,  QY=4.45
 * Extraction                          : QX=4.172, QY=4.23   = wp1

The 2 MHz cavities are for the H0 mode, they are run around 8kV. They are run at the revolution frequency. The beam is below transition.

The 4 MHz cavities are used to flatten the bunches. They reduce the peak current compared to the average currents. The 4 MHz cavities are run between 6-8 kV. They are run at twice the revolution frequency.
