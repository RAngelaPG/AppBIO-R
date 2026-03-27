Version 1.3.14
----------------------------------------------------------------

* Fixed -Wformat, -Wformat-extra-args in Windows oldrelease CRAN checks.

Version 1.3.13
----------------------------------------------------------------

* Added Authors@R field in description to comply with CRAN requirements

* Update to fix a CRAN warning (using non partable flags during compilation)

* Maintainer change. Frederic Bertrand is now the maintainer of the package.
  (https://github.com/wch/Rttf2pt1/issues/25#issuecomment-3320579566).
  
* Fixed C code compilation warnings : -Wcomment, -Wparentheses, 
  -Wunused-function, -Wunused-label, -Wunused-but-set-variable, 
  -Wdeprecated-declarations, -Wuninitialized

* Fixed gcc-ASAN -Wstringop-overflow issue

* Fixed include of <netinet/in.h> on Windows systems

* Fixed -Wimplicit-function-declaration in Windows CRAN checks

* Fixed FIXED issue, -Wmaybe-uninitialized, -Wunused-but-set-variable, 
  -Wunused-function, -Wcomment, -Wparentheses, -Wformat, -Wformat-extra-args 
  in Windows CRAN checks

* Added lib to windows makefile

* Added basic unit tests

Version 1.3.12
--------------

* Closed #10, #19: Reverted some changes from #13 which could cause `ttf2pt1` to crashe, and cause the message `No FontName. Skipping.` to appear. (#20)

Version 1.3.11
--------------

* Closed #15, #17: fixed `-Wformat` and `-Wformat-pedantic` warnings, which raise significant warnings on R-devel. (#16, #18)

Version 1.3.10
--------------

* Fixed signed/unsigned type mismatches.

Version 1.3.9
-------------

* Small fix to stop a compiler warning with gcc11, which raised a WARNING.

Version 1.3.8
-------------

* Small fix to stop a compiler warning with clang 10, which raised a WARNING.

Version 1.3.7
-------------

* Small fix to stop a compiler warning from the `-Wparentheses` flag, which raised a WARNING.

Version 1.3.6
-------------

* Small fix for compilation on Windows.

Version 1.3.5
-------------

* Fixed compilation issues on Windows when it does not have short filename support (at the suggestion of Tomas Kalibera).

Version 1.3.4
-------------

* Fixed compilation on Windows for R 3.3.

Version 1.3.1
-------------

* Made compiler warnings go away, because new version of R CMD check gives warning on compiler warnings.

Version 1.3
-----------

* Changed license line in DESCRIPTION to LICENSE, and addded LICENSE file.

Version 1.2
-----------

* Changed license line in DESCRIPTION to BSD_3_clause.

Version 1.1
-----------

* Force building a 32-bit binary for Windows, by adding the `-m32` flag to gcc in src/Makefile.win, because CRAN builds a 64-bit binary by default, which doesn't work on all systems.

* Added a fallback check for the exec/i386/ directory on Mac, because CRAN puts it there by default.

* Changed the format of the author list in DESCRIPTION.
