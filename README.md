# CS6260-proj
Implementation of Le and Hodgin's Realtime Skinning with Optimized Centers of Rotation in Maya 2017

This project consists of a Maya Deformer plugin of Le and Hodgin's
algorithm for Realtime Skinning using Optimized Centers of Rotation
for Maya 2017. The deformer is in what I would term as an alpha state
and is essentially a prototype for further development. In its current
implementation it consists of a serial CPU and a GPU implementation 
of the deformation algorithm, as well as serial implementation of the
precomputation algorithm, in its naive form.

Reference:
https://www.disneyresearch.com/publication/skinning-with-optimized-cors/

Le, B. H., & Hodgins, J. K. (2016). Real-time skeletal skinning with 
optimized centers of rotation. ACM Transactions on Graphics (TOG), 35(4), 37.

TODO:
Implementation of precomputation using optimizations specified in
Le and Hodgin's paper.

Implementation of a Qt UI to itegrate the deformer into the Maya interface.

Usage:
The plugin has been tested and compiles under Microsoft Visual C++ 2012, 
update 4, on Microsoft Windows operating systems. No testing has been 
performed on Linux or OSX, however, no libraries external to OpenMaya have 
been used, so given proper configurations for projects, the plugin should 
compile on those platforms without issue.

To apply the deformer, select the bones, then mesh and run the bind script.
After weighting your mesh, uncheck valid precomputation and then playback to
force precomputation evaluation.

Disclaimer:
This plugin is not production ready at this time. It is a prototype.  
No warranty is expressed or implied, use at your own risk.

Licensing:
This code is licensce under the Createive Common's Attribution-Non-Commercial-
ShareAlike 4.0 license, details of which can be found here:

https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode

Basically, don't sell it, tweak it to your heart's content and acknowledge 
the use of what you found here. Thanks for playing fair.
