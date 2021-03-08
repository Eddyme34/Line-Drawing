# Line-Drawing
Uses ImGui to draw lines from an input file

Algorithms Implemented
Digital Differential Analyzer(DDA) Line Drawing(line 104) - Uses the algorithm from the textbook with little changes. The line is drawn by incrementing x or y based on dx and dy.
Bresenham Line Drawing(line 138) - Based on the textbook code. It is altered to include other cases like drawing lines to the left, negative cases, and slope greater than 1.
  
  

Translation(line 233) - It creates a 3x3 matrix based on a translation vector input then calls the transformation function which multiplies the matrix with the desired point.


Scale(line 239) - It takes in a point, scale, and center points and used them to create a 3x3 matrix. Then it translated to the origin with the translation function, transforms them, then translates back


Scale(line 252) - It takes in a point, angle, and center points and used them to create a 3x3 matrix. Then it is translated to the origin with the translation function, transforms them, then translates back
  

Line Clipping Cohen-Sutherland and Sutherland Hodgeman(line 320 and 484) - They both use code from the textbook but altered for this program. Cohen-Sutherland’s only difference is that it returns a vector for draw call function and Sutherland-Hodgeman is mostly the same but with bug fixes.
  

Scan-Line Polygon Fill(line 552): This took the longest to implement because I went through several versions of ideas on how to approach this. Some worked(edge and active tables) but crashes the program due to how inefficient it is. But my working approach is iterating through the bounding box for each polygon and changes parity from outside to inside if it intersects an edge. An edge is identified from a boolean buffer created during line drawing. Cases like vertices are solved by creating an edge table and using it to determine if the vertices are connected to y-minimum edges or y-maximum edges and horizontal cases are solved by not including them in the buffer. It’s not perfect however as some angles can trick the algorithm


How to Use: On Windows
I ran the code by compiling it in a CMake Application using Visual Studio Compiler to create the CMake Text file. Then created a launch file to input the text file on Visual Studio Code.
In the GUI for every polygon in an input file, a new tab is created. Each tab represents a polygon you can control separately. The first options are the choice between DDA and Bresenham. The first 4 sliders represent translation vectors, scaling vectors, and angles that move the polygon in real-time. The last 4 sliders are to control the viewport that can clip the polygons. The radio buttons after that are options to rasterize.


Writing to Input File
The button on the bottom of the GUI called “update input” writes all current polygon positions to the input file.
