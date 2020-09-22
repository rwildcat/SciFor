---
project: SciFor  
summary: SciFor -- a modern API to specialized classic libraries   
project_github: https://github.com/rwildcat/SciFor
project_website: https://github.com  
author: Ramon Solano  
email: ramon.solano@gmail.com  
github: https://github.com/rwildcat  
src_dir: ./src  
output_dir: ./doc  
display: public  
         protected  
         private  
source: true  
graph: false  
coloured_edges: true  
exclude: fftw3.f03  
exclude: opkdmain.f  
exclude: opkda1.f  
exclude: opkda2.f  
exclude: src/dev-tests  
predocmark:   
predocmark_alt: >  
md_extensions: markdown.extensions.smarty  
---

Main rutines
-------------

* General      -- Pretty matrices printing
* NumFor       -- Numeric linear spaces, array flipping
* [LAPACK][]   -- Solve `Ax=b`
* [FGSL][]     -- Interpolation
* [FFTW][]     -- FFT and IFFT
* [NetCDF][]   -- Read data, attributes
* [ODEPACK][]  -- Solve sysems of ODEs (`LSODAR()`)
* [RANDLIB90][]  -- Random numbers

To compile
----------

```sh
   $ cd SciFor
   $ cd build
   $ cmake ..
   $ make
   $ sudo make install
```

Documentation
-------------

Internal documentation is included inside the source code, and can be processed using the Fortran Documentator ([ford][]).

To generate:

```sh
$ cd SciFor
$ ford scifor.md
```

Documentation will be generated in the `./doc` directory, as `HTML` pages. Look for the `index.html` file.

### NOTES

   * If `ford` is installed using `pip`, it may lie under the current `python pip bin` directory. For example, on a Mac using `MacPort` to install `Python`:

      `$ sudo pip install ford`

      , it was installed on:
   
      `/opt/local/Library/Frameworks/Python.framework/Versions/3.7/bin/ford`

      To enable `ford` from the commando line, one can just:

      `$ sudo ln -s /opt/local/Library/Frameworks/Python.framework/Versions/3.7/bin/ford /opt/local/bin`

[LAPACK]: http://www.netlib.org/lapack/
[FGSL]: https://doku.lrz.de/display/PUBLIC/FGSL+-+A+Fortran+interface+to+the+GNU+Scientific+Library
[FFTW]: http://www.fftw.org
[NetCDF]: https://www.unidata.ucar.edu/software/netcdf/
[ODEPACK]: https://computing.llnl.gov/casc/odepack/
[Ford]: https://github.com/Fortran-FOSS-Programmers/ford/wiki
[RANDLIB90]: https://biostatistics.mdanderson.org/SoftwareDownload/SingleSoftware/Index/27
[ford]: https://github.com/cmacmackin/ford