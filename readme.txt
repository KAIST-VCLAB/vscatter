--------------------------------------------------------------------------------------------------
README: Multiple Scattering in V-groove - a BSDF plugin of multiple scattering in a V-groove
--------------------------------------------------------------------------------------------------

Multiple scattering in symmetric V-grooves (This version includes multiple importance sampling.)

Written by Joo Ho Lee and Min H. Kim at KAIST VCLAB

Licensed under the GNU GPL license

--------------------------------------------------------------------------------------------------
Background
--------------------------------------------------------------------------------------------------

We provide our implementation of "Practical Multiple Scattering for Rough surfaces", published in SIGGRAPH ASIA 2018.

In this work, we extend the Cook-Torrance BRDF model with multiple scattering, allowing for multiple reflection in a V-groove.

For technical details, refer to our paper:

- Practical Multiple Scattering for Rough Surfaces
  Joo Ho Lee, Adrian Jarabo, Daniel S. Jeon, Diego Gutierrez, Min H. Kim, 
  ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2018)
  http://vclab.kaist.ac.kr/siggraphasia2018p3/index.html
  https://dl.acm.org/citation.cfm?id=3275016

Please feel free to use the following bibtex entry to cite this paper for your academic publication.

@Article{V-Scattering:SIGA:2018,
  author  = {Joo Ho Lee and Adrian Jarabo and Daniel S. Jeon and Diego Gutierrez and Min H. Kim},
  title   = {Practical Multiple Scattering for Rough Surfaces},
  journal = {ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2018)},
  year    = {2018},
  volume  = {36},
  number  = {6},
  pages   = {275:1--12},
  doi     = "10.1145/3272127.3275016",
  url     = "http://dx.doi.org/10.1145/3272127.3275016",
} 

--------------------------------------------------------------------------------------------------
How to compile
--------------------------------------------------------------------------------------------------

Our codes are based on the Mitsuba renderer of version 0.5.0 (https://www.mitsuba-renderer.org/).

Our implementation consists of three files: vscatter.h, roughconductorVgroove, and roughconductorMS.

Before compiling, copy them to your path of the Mitsuba renderer, e.g., mitsuba/src/bsdfs. 

Then modify the "SConsript" file to make scons compile our codes.

Add the following statements in the mitsuba/src/bsdfs/sconscript file:

plugins += env.SharedLibrary('roughconductorMS', ['roughconductorMS.cpp'])
plugins += env.SharedLibrary('roughconductorVgroove', ['roughconductorVgroove.cpp'])

For the test plugin, copy test_ourtest.cpp in mitsuba/src/tests. 

Here, we do not need to modify the sconsript file because scons automatically detects our codes and compiled them already.

Also, test_bsdfours.xml need to be located in the path: dist/data/tests in order to run our test correctly.

--------------------------------------------------------------------------------------------------
About the parameters in our model
--------------------------------------------------------------------------------------------------

Our BSDF plugin has six parameters as follows:

-material: A name of a material preset
-distribution: the type of microfacet normal distribution
-alphaU, alphaV: the anisotropic roughness values along the tangent and bitangent direction
-ScatteringOrderMin: the minimum of reflections in V-groove.
-scatteringOrderMax: The maximum of reflections in V-groove.
-scatteringOrderMaxForPDF: The maximum of reflections for PDF computation.

Note that our model counts k bounces in V-groove (ScatteringOrderMin <= k <= ScatteringOrderMax).

--------------------------------------------------------------------------------------------------
Source Code Files
--------------------------------------------------------------------------------------------------

code/vscatter.h			             : the geometric attenuation functions
code/roughConductorMS.cpp	       : our multiple scattering model with symmetric V-groove
code/roughConductorVgroove.cpp   : our single scattering model with a V-groove G-term.
code/test_ourtest.cpp		         : a test software for our model
data/test_bsdfours.xml		       : the input file of our test evaluation

--------------------------------------------------------------------------------------------------
How to use
--------------------------------------------------------------------------------------------------

1. For quantitative evaluation of energy conservation, reciprocity and single scattering comparison

After compiling the mitsuba renderer, a folder, named as dist, will be created.

To excute the evaluation software, type the following statement in your terminal in the mitsuba folder:

$./dist/mtsutil.exe test_ourtest

2. For qualitative evaluation of multiple importance sampling

Read scene/MIS_test.xml in the mitsuba renderer.

3. For qualitative evaluation of our model (Figure 6)

Read scene/Figure6.xml in the mitsuba renderer.

--------------------------------------------------------------------------------------------------
Version History
--------------------------------------------------------------------------------------------------

 * Version 1.0 (2019-05-20): Initial Release

