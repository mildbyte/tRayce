tRayce
======

Introduction
------------

tRayce is a basic raytracer written in C++, inspired by various articles on the
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
* **Multisample anti-aliasing** (`Scene::doAA`): casts several rays per pixel
    (`Scene::msaaSamples` squared) to make the result smoother. Since using this
    incurs a large performance penalty, there is an option `Scene::msaaOptimize`
    that only performs multisampling at edges of primitives.
* **Soft shadows** (`Scene::softShadowSamples`): when using area lights, instead
    of placing only one point light, approximates the area light with several
    point lights, each contributing a fraction of the total shading, resulting
    in soft shadows.
* **Simple post-processing**: `Scene::render` takes a function that is applied
    to every pixel in the resultant bitmap. It is given the colour of the pixel
    (unnormalised, 0..infinity) and the depth, which allows for some
    effects such as light falloff or playing around with colour intensities.
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

TODO
----

* **Implement irradiance caching/precomputation**: since irradiance values change much slower
    than global illumination, we could cache it at certain points and interpolate during the
    final gather step. Or make every nth photon a radiance photon and query radiance photons
    instead.
* **Stratify**: instead of purely random sampling during final gathering/diffuse reflections,
    use a grid and randomly choose a direction inside the grid, which would eliminate variance.
* **Multicore**: OpenMP?
