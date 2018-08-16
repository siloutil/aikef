/*! @file
 * Defines version number of aikef and contains a small changelog.
 * 
 * The changelog and other comments should be written in a style
 * compatible to <a href="http://www.stack.nl/~dimitri/doxygen/">
 * doxygen </a>.
 *
 * @par Changelog
	 -# Changed Block_Nr from int64 to int32
	 -# Write chosen Statefile type to log  
 * - 6.0: HDF5 Statefile support added & several major bugfixes.
 *   -# Added support for reading and writing HDF5 compatible Statefiles.
 *   -# Removed lots of unused or deprected variables
 *   -# Several major Bugfixes
 *   -# Reworked makefile with support for cluster profiles, parallel compilation and the possibility to start the simulation automatically after compilation
 * - 5.0: Master version of aikef created.
 *   -# Merged the sources of Hendrik Kriegel and Christoph Koenders.
 *   -# Changed data output to HDF5 (integrated compression).
 *   -# Cleaned source code, removed superflous parts, improved
 *      readability and renamed ambigious parts.
 *   -# Created central git-repository for future development at
        http://aikef@aikef.no-ip.org/git/aikef.
 *   -# Tested several compilers and flags for small servers.
 *   -# Changed grid output from `float` to `double`.
 */

#ifndef VERSION_H
#define VERSION_H


/*! Major version number of aikef.
 *
 * This version number should be incremented if major changes in the source occur
 * which might break backward compatibility or new features are added.
 */
#define AIKEF_VERSION_MAJOR 6




#endif /* VERSION_H */

