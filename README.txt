README.txt, last updated August 9, 2011 (Karl Yando)

Description
-----------

MExED_v1.4 is a Geant4-based simulation of the second generation MEPED electron telescopes flown aboard the NOAA/TIROS-N (POES) weather satellites.  This folder contains the C++ code for the model geometry and simulation physics, and may be compiled with a working Geant4 installation.  


Known Issues and Assumptions
-----------

1) During program exit, a segmentation fault occurs, outputting support information to standard error.  This error occurs only after simulation is complete, and does not appear to affect the quality of simulation output.  
2) Standard Geant4 volume intersection analysis was performed using the recursive grid test, and the following volumes were tweaked:
- Contacts 2 and 3: No overlap should exist here, but the volumes were foreshortened so as to leave a thin (2 micron) gap.
3) The following simplifications were made with respect to materials:
- the detector ring mounts were modeled as Kel-F insulator instead of glass expoxy.


Version History
-----------
v1.4            Implemented April 2010.  Revised output format.  Public release August 2011.
        Class "G4MultipleScattering" replaced with particle specific classes per community recommendation.
v1.3            Implemented October 2008.  Standardization of codebase and fork of MExED.  
v1.2            Implemented May 2008.  Basis for undergraduate thesis.
pre-1.2 Developmental use only.


Authors
-----------
Developed at Dartmouth College from example "N03" by Karl Yando (2007-2011), under the direction of Professor Robyn Millan.  

The simulation physics employed in this program were written by Dennis H. Wright (SLAC), and may be obtained in unmodified form from the following URL:
http://www.slac.stanford.edu/comp/physics/geant4/slac_physics_lists/micro/physlistdoc.html

The simulation physics files are:
SEBosonPhysics.cc
SEDecayPhysics.cc
SEHadronPhysics.cc
SEIonPhysics.cc
SELeptonPhysics.cc
SENeutronPhysics.cc
SEPhysicsList.cc
(and their corresponding header files)


Contacts
-----------
please report issues to the appropriate GitHub repository


Acknowledgements
-----------
This project made possible by:

    The Geant4 Collaboration (http://geant4.cern.ch/).

    Dennis H. Wright (SLAC), author of the Geant4 physics list used herein.


License
-----------
This compilation is released under an "Open Research License" (ORL), written in the spirit of the ORL proposed by Victoria Stodden.  Various components are released under either the 2-clause (Simplified) BSD license or the Creative Commons Attribution (BY) 3.0 license.  Any code derived from the Geant4 toolkit is subject to the Geant4 Software License, and copyrights for any third-party code rest with their respective rightsholders. 

Code components are served by the 2-clause BSD license:
-----------
Copyright (c) 2007-2011, Karl Yando.  All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
-----------
Text, documentation, figures, and other media are served by the Creative Commons Attribution 3.0  license (CC BY 3.0), a copy of which may be found at the URL:
http://creativecommons.org/licenses/by/3.0/
-----------

