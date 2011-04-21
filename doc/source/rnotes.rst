Release Notes
=============

Version 0.4.x
-------------
 * Added new robust statistical methods
 * Added Global Sky Model models to skymap
 * Added new smoothing and fitting fuctions to misc.mathutil
 * Added a progress bar for use with long-running applications
 * Added in new LWA-1 dipole NEC models
 * Re-worked the plotAntenna.py script to use the pre-computed dipole response fits
 * Added support for EXCITATION keywords in misc.nec_util
 * Added new gofast readers for TBW, TBN, and DRX
 * Added better polyphase filter support
 * Fixed initiaization problems the various lsl.reader Frame objects
 * Fixed the TBN ring buffer
 * Added support for simulating point sources in TBW and TBN data
 * Moved the tests directory out of lsl proper
 * Added new offset/average options to drxSpectra and tbnSpectra
 * Added new DRSU direct access module for linux system (experimental)
 * Added documentation for the various C extensions

Version 0.3.x
-------------
 * New setuptools-based build system
 * New, improved documentation
 * Added post-acquisition beam forming module
 * Added in modules for writing simulated data to the S60, TBW, TBN, and 
 * DRX raw data formats
 * Added a module for simulated basic test signals in TBW, TBN, and DRX data
 * Moved the lwa_user script 'astrostatus.py' into the scripts directory
 * Fixed the FITS IDI writer
 * Fixed a couple of bugs dealing with incorrect chunking of data in scripts/tbwSpectra.py and scripts/tbnSpectra.py

Version 0.2.x
-------------
 * Added in other modules from the lwa_user library developed by NRL

Version 0.1.x
-------------
 * Initial version
