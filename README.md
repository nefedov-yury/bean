---------------------------------------------------------------------------
## IMPORTANT

Before running cmake

0. Bean requires at least [CMake](https://cmake.org/download) 3.15

1. [CERN ROOT](https://root.cern.ch) must be installed.\
   The `root-config` script must be in `$ROOTSYS/bin` directory,
   otherwise, you must specify the path to `root-config` when call
   the cmake.
   For example, in IHEP cluster (CentOS7):
   ```
   cmake -DROOT_CONFIG_SEARCHPATH="/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/lcg/LCG_84/ROOT/6.20.02/x86_64-centos7-gcc49-opt/bin"
   ```

2. [CLHEP](https://proj-clhep.web.cern.ch/proj-clhep) library must be
   installed.\
   The path to CLHEP can be specified in the environment variable
   `CLHEP_DIR`. Also you can specify the path to CLHEP when calling
   cmake.
   For example, in IHEP cluster (CentOS7):
   ```
   cmake -DCLHEP_SEARCHPATH="/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/lcg/LCG_84/clhep/2.4.4.0/x86_64-centos7-gcc49-opt"
   ```

3. Some notes for installing the BEAN on Windows are
   [here](doc/Win_readme.md)

---------------------------------------------------------------------------
## INSTALLATION

1. Create build directory: `mkdir build; cd build`

2. Run `cmake ../ [cmake-options]` to prepare the infrastructure
   for building the BEAN program.
   * For example:\
     `cmake ../ -DBOSS_VERSION=6.6.4`\
     or\
     `cmake ../ -DBOSS_VERSION=7.0.4`
   * If you want to compile program using a specific compiler,
     use the following cmake options:\
     `-DCMAKE_C_COMPILER=gcc49 -DCMAKE_CXX_COMPILER=g++49`
   * More detailed printing can be achieved by adding the "log-level"
     option to cmake:\
     `cmake --log-level=DEBUG`
   * If you wish, you can run `cmake-gui` and change the options
     in the dialog interface.

3. Compile the BEAN by runing `make` and then `make install`
   * To see detailed output of compilation commands, use
     `make VERBOSE=1`\
     Alternatively, when calling cmake you can use the option:
     `-DCMAKE_VERBOSE_MAKEFILE=ON`

   * `make install` will install the `bean_version.exe` into
     the `workdir/` directory, and the libraries into the
     `workdir/BeanLib_version` directory. Here the version is the
     number you specify as BOSS_VERSION when you call cmake.

4. The BEAN program uses sqlite database located in the
   `Analysis/DatabaseSvc/dat/` directory.
   **Currently, updating this database is only possible "manually".**
   You need to go to the `https://docbes3.ihep.ac.cn/db_ana/` website,
   download and unzip two files: `offlinedb.db` and `run.db`.
   + _Note: You can check the integrity of downloaded files using
     the `db.timestamp` file, which contains the Unix time and
     md5 checksums for those files._

5. To remove BEAN files from the `workdir/` directory, use the
   `make uninstall` command. This may be useful when updating
   the program. Alternatively, you can delete the 'bean_version.exe'
   file and the 'BeanLib_version' directory manually.

6. After running make install, you can remove the build directory
   at any time.

---------------------------------------------------------------------------
## NOTES AND NEWS

* The user source programs MUST be in directory `BeanUser/`.\
  Edit `BeanUser/CMakeLists.txt` to add name of your file
  to the list of files to be compiled.
  On Windows, be sure to include the `DLLDefines.h` header and use
  the `BeanUserShared_EXPORT` macro before your functions.

* The AbsCor algorithm in BEAN does not use the database, unlike
  the BOSS version. However, you can read calibration files using
  the functions:
  ```
    void ReadDatac3p(std::string DataPathc3p, bool info=true);
    void ReadParMcCor(std::string paraPath, bool info=true);
    void ReadCorFunpara(std::string CorFunparaPath, bool info=true);
  ```
  The path to the calibration files in CernVM `/cvmfs/` file system
  can be obtained using the python3 program `GetAbsCorFiles.py` located
  in `bean/scripts` directory.\
  An example of using these functions is contained in the file:
  `BeanUser/TestAbsCor.cxx`

* A small bug in the VertexFit package where the total energy was not
  updated after fitting has been fixed for all versions of the BEAN
  program.

* The `DecayTable` class has been added to `EventTag` algorithm.
  This class is designed to fill out a Monte Carlo particle decay
  table in text form.

* The `RscanDQ` package has been added to the BEAN algorithms.
  When you using the 2015 Rscan data, you must use this package
  to exclude bad runs. For more information, see the
  [Tau and QCD Group Data Samples page](
  https://docbes3.ihep.ac.cn/~tauqcdgroup/index.php/Data_Samples).

* The `TrackCorrection` class for the helix parameters corrections
  for MC tracks, based on "TrackCorrection" package has been added
  to the BEAN algorithms.
  For more information, see the [Charmonium Group's Data Quality page](
  https://docbes3.ihep.ac.cn/~charmoniumgroup/index.php/DataQuality_Page).

---------------------------------------------------------------------------
## DOCUMENTATION

* Examples in the `BeanUser/` directory:

| File             | Description                                       |
| :---             | :---                                              |
| ---------------- | ------------------------------------------------- |
| UserTest.cxx     | two simple examples of writing user functions in  |
| User1.cxx        | the BEAN; demonstration of saving histograms      |
| ---------------- | ------------------------------------------------- |
| Rhopi.cxx        | Rhopi program from `Analysis/Physics/RhopiAlg`    |
| ---------------- | ------------------------------------------------- |
| TestAbsCor.cxx   | example of using functions of AbsCor algorithms   |
| TestPID.cxx      | example of using the ParticleID library functions |
| MagField.cxx     | example of using the MagneticField functions      |
| TestEventTag.cxx | example of using the EventTag facility            |
| TestDb.cxx       | example of using the SqliteInterface class        |
| ---------------- | ------------------------------------------------- |
| JobInfo.cxx      | program for printing information from JobInfoTree |
| EntryList.cxx    | program with an example of using TEntryList       |
| ---------------- | ------------------------------------------------- |

* Note that the BEAN documentation on the
  [BES3 Offline Software Group page](
  https://docbes3.ihep.ac.cn/~offlinesoftware/index.php/BEAN)
  is outdated and does not seem to be very useful.

---------------------------------------------------------------------------
### ROOT-PROOF: Deprecated and no longer supported

* The PROOF framework is marked as deprecated in the ROOT since
  version 6.26, and was completely removed in version 6.32.
  Therefore, the BEAN was rewritten so that it can be compiled
  without the PROOF support.

* If you need PROOF functionality, the old 'CMakeLists.txt' and
  other PROOF related files are in 'doc/OldStuff/PROOF/' directory.
  BEAN running in PROOF mode has not been tested for a long time,
  so use it at your own risk.

---------------------------------------------------------------------------
### TODO and legacy notes

* Add an example of using the `DecayTable` class in `TestEventTag.cxx`

* `make updatedb` command was intended to be used to install or update
  a database. After entering password access, this method no longer
  (or not yet?) works.

* `make setup_file` command creates the `setup_BOSS_VERSION.sh` file.
  This bash script sets environment variables with paths to libraries used
  in the BEAN. In the case of a standard installation of all components,
  this script **IS NOT REQUIRED**.
  Perhaps it could be useful for debugging.

---------------------------------------------------------------------------
