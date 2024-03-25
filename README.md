# INF585 - Computer Animation - Final Project - Creation of a 3D physics simulation scene

## Description

This project was realized as part of the course INF585 - Computer Animation at Ecole Polytechnique.
In this project, I tried to implemented the sand modelisation paper _Animating Sand as a Fluid_ by Yongning Zhu and Robert Bridson in the ACM Transactions on Graphics 24(3):965-972 of 2005.

## Installation

To install and run, just clone the archive and then follow the compiling steps relative to your OS given on this [page](https://imagecomputing.net/cgp/compilation/content/01_compilation/index.html).

## Controls in the project

Once the project is launched, the mouse commands are displayed in the console.
I also display them here:

-  Mouse left click + drag: Rotate the camera (pitch/yaw) around its focus point.
-  Mouse right click + drag: Camera move close/far from the central focus point (the focus point remains unchanged).
-  Ctrl + Mouse left click + drag: Translate/Pan the camera and its central focus point in the viewspace plane.
-  Ctrl + Mouse right click + drag: Translate the camera and its central focus point in front/back direction.
-  Shift + Key left/right (or key r/f): Rotate the "up" direction used in this Euler angle representation (rotation around z).

The other commands are available inside the User Interface that is displayed on the visualization (you might have to click on it to open it first).
It is possible to toggle what is displayed (the actual particles, the bigger particles that make the rendering better, the grid (only the outside lines, to give an idea of the size of the cells), the reconstructed surface and the walls). It is also possible to activate or deactivate the correction of the velocity using the divergence, to modify the position of the bottom-left angle of the pile of sand, to restart the animation and to modify the time steps between frames.

## Results

Here are some key frames displaying surface reconstruction:

![](/media/start_anim.png)

![](/media/middle_anim.png)

![](/media/end_anim.png)

Here is what the display looks like when using the bigger radius for better result:

![](/media/with_part.png)

There are some incoherencies in the movement of the sand that I was not able to solve:

![](/media/layers.png)

And you can download the video of the final animation [here](/media/project_video.mp4).

The report for the project is available [here](/media/final_report.pdf).

## Author and aknowledgment

Author: Marie Audouard, using the [CGP library](https://imagecomputing.net/cgp/01_general/index.html) created by Damien Rohmer.
