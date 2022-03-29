# PS Booster Optics Repository
Optics repository for PS Booster Accelerator

## layoutDB directory

Contains the files as downloaded from the Layout Database

Steps to follow for an update to a new layout:
* download the new file from:  https://layout.web.cern.ch/madsequence.aspx?machine=PSB
* Split the file to produce the .seq, .ele, and .str files using the script SplitLayoutDB.py in the cmd directory

> > pyhton ../cmd/SplitLayoutDB.py -i <layoutDB filename> -o psb

* compare the new files with the old ones to identify differences
* copy the new files .seq, .ele and .str to the corresponding directories. Use the scrip in the cmd dirctory

** Layout Changes (recent history) **

* 17.03.2015 : Finemet cavity installed: BR4.ACWF6L1  
* New definition of collimator: ECOLLIMATOR and RCOLLIMATOR are now defined with the help of an APERTURE command. E.g.  COLLIM: ECOLLIMATOR, L = 0.5, APERTYPE = ELLIPSE, APERTURE = 0.01,0.005;
* Changes during EYETS16-17 are documented in: EDMS no 1803330 https://edms.cern.ch/document/1803330/0.01
* Changes during YETS17-18 :
* Changes during LS2:  
   * new injection in period P01 
