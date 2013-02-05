sdl-vertex-animation

Author:

	Jurgen Schober
	
Date:
   
	February, 2013
	
Short:
  
	Example using vertex arrays and vertex animation with OpenGL.

Description:

	Based on sdl-render-to-texture this example uses code from previous sdl-xx-examples I've written.
	This example focuses on using vertex array and vertex manipulation to generate a waving flag.
	
	All sdl-xx-examples are written in C++11 using MinGW gcc 4.6 and are Windows only. I'm using
	Eclipse Juno as Development IDE.

Libs used:

	boost_thread
	boost_system
	boost_filesystem
	boost_iostreams
	glew
	+ OpenGL

License:

	Use as is. No license other then the ones included with third party libraries are required.

Compiler used:

	MinGW with Eclipse Juno. Windows only. Linux might just work, MacOS will need some work due to 
	the fact OSX needs to run the render loop in the main loop. This example runs a render thread decoupled
	from the main thread.

Have fun
Jurgen
