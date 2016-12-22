# AF_randomc_h
Header only version of the part of [Agner Fog](http://www.agner.org/)'s [randomc](http://www.agner.org/random/randomc.zip) library (versioned 2014-Jun-14).
Contains the CRandomMersenne, CRandomMother and CRandomSFMT classes with some small (mostly technical - just to make MSVC2015 happy) corrections. This classes implements very fast random number generators with a good randomness.

Made for the [NNTL](https://github.com/Arech/nntl) project, but distributed as a standalone module due to license incompability
(randomc is licensed under the GNU General Public License http://www.gnu.org/licenses/gpl.html)

If you want to use it outside of the NNTL, you may probably want to download it from the [original](http://www.agner.org/random/) page.

#### Note for overclockers

The CRandomSFMT class contains very tight loops of SSE instructions which may lead to hardware faults / BSODs on overclocked CPUs. I had to lower a CPU frequency from 3900Mhz to 3500Mhz and raise some voltages to make my PC stable again (it had been working fine on 3900Mhz for a half of a year before that, successfully passing all available pretty intensive stability tests). You've been warned.
