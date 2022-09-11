<h1 align="center">
ulmo
</h1>
<h4 align="center">
Model of Ekman heat transport for slab oceans, now in modern Fortran!
</h1>
Instructions
</h1>
->Pull the repository
</h1>
->Make a new python environment (maybe call it ULMO)
</h1>
->Conda install gfortran
</h1>
->To run the code use "Make run"
</h1>
->"Make clean" removes complied files
</h1>
->"Make all" compiles files only
</h1>
->Version propmt after running refers to the heat transport mechanisms enabled
</h1>
-> 1 = no transport
</h1>
-> 2 = Diffusion
</h1>
-> 3 = Diffusion and Ekman transport
</h1>
-> Constants such as the Diffusion constant can be changed in the constants.F90 file
</h1>
-> Examples can befound and run from the examples folder, (these have pre determined versions and are named accordingly)
</h1>
IF USING THE FGSL VERSION ONLY:
</h1>
->Make a new python environment
</h1>
->Conda install fgsl
</h1>
->Edit make file to compile the correct files e.g ulmo_fgsl.F90 and change the compiler flags



</h4>

<p align="center">
  <img src="https://img.shields.io/badge/wip-%20%F0%9F%9A%A7%20under%20construction%20%F0%9F%9A%A7-yellow"
       alt="wip">
</p>

<p align="center">
<a href="https://fortran-lang.org/">
<img src="https://img.shields.io/badge/fortran-2003-purple.svg"
     alt="Fortran 2003"></a>
<a href="https://www.python.org/downloads/">
<img src="https://img.shields.io/badge/python-3.10-blue.svg"
     alt="Python 3.10"></a>
<a href="https://github.com/psf/black">
<img src="https://img.shields.io/badge/code%20style-black-000000.svg"
     alt="black"></a>
<a href="LICENSE">
<img src="https://img.shields.io/badge/license-MIT-green.svg"
     alt="License: MIT"></a>

<a href="https://github.com/exoclim/ulmo/graphs/contributors">
  <img src="https://img.shields.io/github/contributors/exoclim/ulmo"
       alt="Contributors">
</a>
</p>

