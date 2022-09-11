<h1 align="center">
ulmo
</h1>
<h4 align="center">
Model of Ekman heat transport for slab oceans, now in modern Fortran!
</h4>
Instructions <br />  
->Pull the repository <br />
->Make a new python environment (call it "ulmo") <br />   
->Conda install gfortran  <br />
->To run the code use "Make run" <br />  
->"Make clean" removes complied files  <br />
->"Make all" compiles files only  <br />
->Version prompt after running refers to the heat transport mechanisms enabled <br />  
-> 1 = no transport  <br />
-> 2 = Diffusion  <br />
-> 3 = Diffusion and Ekman transport<br />  
-> Constants such as the Diffusion constant can be changed in the constants.F90 file<br />  
-> Examples can be found and run from the examples folder, (these have pre determined versions and are named accordingly) <br />  
IF USING THE FGSL VERSION:  <br />
->Make a new python environment  <br />
->Conda install fgsl and gfortran <br />
->Edit make file to compile the correct files e.g ulmo_fgsl.F90 and change the compiler flags <br />




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

