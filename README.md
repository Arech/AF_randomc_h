# AF_randomc_h
Header only version of a part of [Agner Fog](http://www.agner.org/)'s [randomc](http://www.agner.org/random/randomc.zip) library (versioned 2014-Jun-14).
Contains CRandomMersenne, CRandomMother and CRandomSFMT classes with small (mostly technical - just to make MSVC2015 happy) corrections. This classes implements very fast random number generators with a good randomness.

Made for [NNTL](https://github.com/Arech/nntl) project, but distributed as standalone module due to license incompability
(randomc is licensed under GNU General Public License http://www.gnu.org/licenses/gpl.html)

If you want to use it outside of NNTL, you probably may want to download it from [original](http://www.agner.org/random/) page.

#### Note for overclockers

CRandomSFMT class contains very tight loops of SSE instructions which may lead to hardware faults / BSODs on overclocked CPUs. I had to lower CPU frequency from 3900Mhz to 3500Mhz and raise some voltages to make my PC stable again (it had been working fine on 3900Mhz for half a year before that, successfully passing all available pretty intensive stability tests). You've been warned.
