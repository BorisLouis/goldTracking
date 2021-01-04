# goldTracking
Algorithm to track gold nanoparticle self assembly with known number of particles


# mainTracking
The main code to be used is called mainTracking.m. It can currently track up to 8 particles and is expecting fluorescence/darkfield type data (bright foreground).
It can currently process .mp4 and .ome.tif as well as entire folder of the two (processing every file in the folder. It is expected to be given a folder in which a single file is OR a folder in which subfolder each containing single files are.

To run, the user needs to provide several input:
1) Path to the folder to analyze (UI)
2) delta - radius around the particles to be picked in pixel. This is to reduce the processing time, smaller is better but you need to be sure to not crop smaller than where you particles are moving
3) number of particles - This code is made to analyze and always find particles, but it needs to receive the correct number of particles preamptively
4) width - width to be used for the gaussian fitting 3 for 200 nm beads and 95nm pixel for instance. 0 Will let the width a free parameter
5) pxSize - pixel size in nm
6) minDist - minimum distance expected between particles in pixels
7) scaleBar - size of the scale bar for gif in um
8) tail - number of frame of the tracked trace to be plotted on the gif movie for every single frame
9) frameRate - frameRate of the movie
10) info.type - 'normal' or 'transmission' transmission will make the code revert the colorscale to create a bright foreground image
11) OutputFolder - name of the folder to be created to store the results. The folder will be created in the path of the file currently analyzed
