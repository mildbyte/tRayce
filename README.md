tRayce
======

Introduction
------------

tRayce is a basic raytracer written in C++, inspired by various articles on the
Internet (such as [this one](http://www.codermind.com/articles/Raytracer-in-C++-Introduction-What-is-ray-tracing.html)).

It supports basically nothing, including, but not limited to:

* Not having a scene definition language, unless we call C++ a scene definition language.
* Not supporting primitives other than spheres, planes and axis-aligned boxes.
* Not supporting changing the camera's direction. Meaning axis-aligned boxes
    don't exactly look exciting.
* Not using CUDA/SSE.
* Not being able to export the rendered image into anything but an uncompressed
    BMP file.
* Not supporting textures, procedural or not, unless we call a plain colour a texture.
* Not being realtime, unless we call 30 FPS at 160x100 on a 2.1 GHz dualcore realtime.
* Not using any acceleration structures to make the rendering faster.
* As of August 28th, it has also been certified to be slow as hell on a Raspberry Pi Model B: ~80 times slower for the same 1280x800 scene, probably because of the kernel that U have there not using hard-float ABI.

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
    and `width`), the focus distance (`Scene::camera.planeDistance`) and the
    position (`Scene::camera.position`). The direction of the camera can't be
    adjusted yet.
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
