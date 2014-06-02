tRayce
======

Introduction
------------

tRayce is a basic raytracer, photon mapper and path tracer written in C++, inspired by various articles on the
Internet (such as [this one](http://www.codermind.com/articles/Raytracer-in-C++-Introduction-What-is-ray-tracing.html)).

It supports basically nothing, including, but not limited to:

* Not having a scene definition language, unless we call C++ a scene definition language.
* Not supporting primitives other than spheres, planes and axis-aligned boxes.
* Not using CUDA/SSE.
* Not being able to export the rendered image into anything but an uncompressed
    BMP file.
* Not supporting textures, procedural or not, unless we call a plain colour a texture.
* Not being realtime, unless we call 30 FPS at 160x100 on a 2.1 GHz dualcore realtime.
* Not having a good enough photon mapping (it kind of casts them and stores them, but it only looks nice with a final gathering step with horrendously large amounts of final gather rays, which is rather slow since I don't do any irradiance caching (at the time of writing, it's been rendering a 320x240 scene with 175000 photons, 500 photons used in an irradiance estimate and 512 final gather rays per pixel) for about a couple hours now)
* Not using any acceleration structures to make the rendering faster.
* (that's not exactly true anymore, I think I have a kd-tree going on there for storing photons during photon mapping)
* As of August 28th, it has also been certified to be slow as hell on a Raspberry Pi Model B: ~80 times slower for the same 1280x800 scene, probably because of the kernel that I have there not using hard-float ABI.

Compilation
-----------

    make

I wrote tRayce on a Windows system and then managed to compile it on Linux
without any major changes (except for rewriting the bitmap export code). The
Makefile in the root folder seems to work, though I haven't tested it on anything
but my system.

You do need a compiler that supports the C++11 standard.

Usage
-----

There is no scene definition language as of yet: the scene is defined in the
source code and compiled with the raytracer. The sample scene already in
`main.cpp` essentially showcases everything tRayce can do for now.

Some features that tRayce does support:

* **Adjusting the camera**: `Scene.camera` sets up parameters such the dimensions
    of the image plane the rays would be cast through (`Scene::camera.height` 
    and `width`), the focus distance (`Scene::camera.planeDistance`), the position
    (`Scene::camera.position`) and look-at direction (`Scene::camera.direction`).
* **Simple post-processing**: `Scene::render` takes a function that is applied
    to every pixel in the resultant bitmap. It is given the colour of the pixel
    (unnormalised, 0..infinity) and the depth, which allows for some
    effects such as light falloff or playing around with colour intensities.
* **Multithreading!**
* **Spheres!** **Triangles!** **Planes!** **Axis-aligned bounding boxes!**
* **Specific to raytracing**:
    * **Multisample anti-aliasing** (`Scene::doAA`): casts several rays per pixel
        (`Scene::msaaSamples` squared) to make the result smoother. Since using this
        incurs a large performance penalty, there is an option `Scene::msaaOptimize`
        that only performs multisampling at edges of primitives.
    * **Soft shadows** (`Scene::softShadowSamples`): when using area lights, instead
        of placing only one point light, approximates the area light with several
        point lights, each contributing a fraction of the total shading, resulting
    in soft shadows.
* **Photon mapping**: After me skimming through dozens of various articles and papers,
    it kind of works and looks nice and even uses a kd-tree to store photons and perform
    nearest-neighbour queries in logarithmic time. What bothers me is that it
    doesn't actually look logarithmic, since it appears that the amount of kd-tree nodes
    visited per pixel grows linearly with the number of stored photons. I believe it has
    something to do with inefficient median splits during the construction (since the tree
    is balanced at build-time and the lookup routine has to access both branches of the tree
    fairly often, meaning that the median is closer to the target point than our current best
    and we are not sure there isn't anything better in the other branch). Also, there is
    something wrong with the brightness of the image: using classic raytracing results in
    much brighter images (by a factor of a thousand) than the same images rendered using a
    photon map. Final gathering is quite slow (which is expected) and sometimes has weird
    artifacts (like white pixels that could not have appeared from anywhere).

    Irradiance "caching" is now supported: instead of computing irradiance every time during
    the final gather, the mapper now looks for the nearest irradiance photon (a photon on the
    position of which the radiance estimate has already been performed). This considerably
    speeds things up. It also looks like a cool Voronoi diagram when rendered directly.

    Using Monte Carlo integration with stratified samples and a Mersenne Twister
    RNG, I finally managed to achieve an almost noiseless 320x240 image with 50000 initial
    photons, 500 photons used for the irradiance estimate and 64x64 final gather samples.
    The image was rendered in ~2 hours. That seems much slower than what is claimed in the
    papers (on their rather ancient hardware)

    Halton sequences, which provide a low-discrepancy sequence that sort of looks like
    stratified random sampling, are also implemented. Setting `Scene::samplingMode = HALTON`
    enables them, with `photonGatherSamples` now being the number of elements in the
    sequence that are used (unlike the `STRATIFIED` mode that uses the square of that number
    of samples). With low amounts of samples, it looks quite bad, with blobs of light smeared
    all over the objects, but with ~1000 it begins to fade away and only be noticeable on
    close inspection.

    Photon visualisation is also supported, but that's mostly for confirming they are landing
    in the right areas.

    Photon maps can be saved and loaded, but they are not checked for validity during load time.
* **Pathtracing**:
    Seems to be what I like working on: simple raytracing is just isn't pretty anymore!
    * Enabled by scene.renderingMode = PATHTRACING;
    * `pathTracingSamplesPerPixel` does what it says on the tin. `pathTracingMaxDepth` limits the
        trace depth. Note that when a ray goes through an object, that's 2 "bounces": when it comes
        in and when it comes out.
    * `Camera::focalDistance` determines the position of the focal plane, whereas `Camera::lensRadius`
        determines the radius of the circle from which the ray will be cast through the focal plane.
        This simulates the depth-of-field effect.
    * Can get some noisy caustics

TODO
----

* **Speedups**: irradiance caching; figure out why I need so many final gather samples. Also,
    there is a nice paper called Balancing Considered Harmful that proposes a different way of
    organising the kd-tree to avoid unnecessary visits to the other branch. Another way is
    proposed in the paper called It’s okay to be skinny, if your friends are fat.
* **Improve the pathtracer**: bidirectional path tracing, then probably move on to Metropolis
    Light Transport
* **Split** the pathtracer, photon mapper and the raytracer
* **More speedups**: if I want to render more complex scenes that have tons of triangles,
    a good data structure is vital.

References
----------
* Photon mapping: Global Illumination using Photon Maps by Henrik Wann Jensen, 
    http://graphics.ucsd.edu/~henrik/papers/photon_map/global_illumination_using_photon_maps_egwr96.pdf 
    The 2000 SIGGRAPH course, "A Practical Guide to Global Illumination using Photon Maps",
    http://graphics.stanford.edu/courses/cs348b-00/course8.pdf
    and lecture notes on Realistic Image Synhesis: Photon Mapping by Philipp Slusallek,
    Karol Myszkowski and Vincent Pegoraro, http://graphics.cs.uni-saarland.de/fileadmin/cguds/courses/ss10/ris/slides/RIS11PhotonMap.pdf
* Irradiance photons: Faster Photon Map Global Illumination by Per H. Christensen,
    http://www.seanet.com/~myandper/jgt99.pdf. Claimed "Rendering this
    image at resolution 1024 × 1024 pixels with up to 16 samples per pixel took 5 min-
    utes and 29 seconds. Out of this time, computing the soft indirect illumination
    took 2 minutes and 45 seconds. There were 11,800 final gathers with a total of
    6.5 million final gather rays." on a "233 MHz Pentium processor and 32 megabytes of memory",
    so my implementation is really slow.
* KD-trees: An introductory tutorial on KD-trees by Andrew W. Moore,
    http://www.autonlab.org/autonweb/14665/version/2/part/5/data/moore-tutorial.pdf?branch=main&language=en
* !A message on the Blender mailing list about trying to implement final gather and having the
    same problems as me: http://lists.blender.org/pipermail/bf-committers/2009-April/023044.html
