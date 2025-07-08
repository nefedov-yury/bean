---------------------------------------------------------------------------
## My experience installing BEAN on Windows
* Tested on Windows 10 version 22H2
* The [Windows Package Manager command-line tool](
  https://learn.microsoft.com/en-us/windows/package-manager/winget/)
  `winget` is used to install a number of programs

---------------------------------------------------------------------------
### For ease of use

1. Install and get started setting up [Windows Terminal](
   https://learn.microsoft.com/en-us/windows/terminal/install)

   `PowerShell> winget install -e --id Microsoft.WindowsTerminal`

2. Git Bash

   It's part of the [Git For Windows](https://gitforwindows.org/) package:
   "The BASH emulation behaves just like the `git` command in LINUX and
   UNIX environments"

   `PowerShell> winget install -e --id Git.Git`

   * Enable git-bash integration with Windows Terminal

---------------------------------------------------------------------------
###  Prerequisites

3. [CMake](https://cmake.org/download)
   * The build has been tested with cmake-3.29, but should work with
     versions higher than 3.15

   `PowerShell> winget install -e --id Kitware.CMake`

4. [Visual Studio Community 2022](
   https://visualstudio.microsoft.com/vs/community)

   `PowerShell> winget install -e --id Microsoft.VisualStudio.2022.Community`

   * Eneble integration with Windows Terminal: ‘Developer PowerShell
     for VS 2022’ or ‘Developer Command Prompt for VS 2022’

   * It is more convenient to use `git bash` with Visual Studio
     compiler in Windows Terminal.
     To do this, edit the Windows Terminal configuration file
     (open JSON file in the parameters of Windows Terminal)
     and add the following lines:
     ```
     {
        "commandline": "\"C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat\" && \"C:\\Program Files\\Git\\bin\\bash.exe\" --login",
        "guid": "{c4a59547-3fc1-4ba9-a0c0-5a81bcc7f9e6}",
        "hidden": false,
        "name": "Git Bash with VS 2022 CLI"
     },
     ```
     + _NOTE1: change the program paths according to your settings_
     + _NOTE2: a new guid can be generated with:_
       `PowerShell> [guid]::NewGuid()`

5. [Root](https://root.cern.ch)

   * root version 6 is preferred, we used version root-6.30
   * download the binary built with the exact same version
     of Visual Studio than the one installed on your system
   * use exe-installer: a regular Windows installer package,
     that also sets up the required environment variables
   * check that the `ROOTSYS` and `PATH` environment
     variables are correct:\
     `PowerShell> echo $Env:ROOTSYS`\
     `PowerShell> echo $Env:PATH`

6. [CLHEP](https://proj-clhep.web.cern.ch/proj-clhep/clhep23.html)
   or [CLHEP git repository](https://gitlab.cern.ch/CLHEP/CLHEP)

   * The latest available version is preferred, we used clhep-2.4.7.1

   * Unfortunately, building for the modern version of the
     Visual Studio C++ compiler is not supported.
     Therefore, a number of changes need to be made:
     - In file `cmake\Modules\ClhepVariables.cmake`
       change flags for C++ compiler:
       ```
       set(CMAKE_CXX_FLAGS "/nologo /Zc:__cplusplus /std:c++17 /MD /GR /EHsc /D USING_VISUAL")
       ```

     - In a number of files in the `Random\test\` directory,
       it is necessary to replace the gnu-specific attribute:
       `__attribute__ ((unused))` with a standard C++17 attribute:
       `[[maybe_unused]]`

       For example, the following changes were made to the file
       `Random\test\testAnonymousEngineRestore.cc`:
       ```
       // double __attribute__ ((unused)) r = 0;
       [[maybe_unused]] double r = 0; // C++17
       ```

   * The following commands are used to build CLHEP:
     ```
     GitBashVScli> cmake ../ \
     -G "NMake Makefiles" \
     -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
     -DCMAKE_INSTALL_PREFIX="${INSTDIR}" \
     -DCLHEP_BUILD_DOCS=OFF
     GitBashVScli> nmake
     GitBashVScli> nmake install
     ```
     + _NOTE: change the path to the directory where you want
       to install CLHEP, for example:_
       `INSTDIR="${HOME}/CLHEP/win64_vc17/"`

   * Please define an environment variable `CLHEP_DIR` containing
     the path to the installed version of CLHEP.
     ```
     PowerShell> [Environment]::SetEnvironmentVariable( "CLHEP_DIR", 'C:\Users\<USERNAME>\CLHEP\win64_vc17', 'User')
     ```
     This variable is used in cmake to build the BEAN program.

---------------------------------------------------------------------------
### Build BEAN

* To build the program, a standard set of commands is used:
  `cmake`, `nmake`, `nmake install`

  _NOTE: The `PreLoad.cmake` file specifies that the "NMake Makefiles"
  generator should be used on Windows._\
  **The use of other generators, for example Visual Studio 17 2022,
  has not been tested and most likely will not work.**

---------------------------------------------------------------------------
