// Define characteristic length
lc = 0.5;

// Define points
Point(1) = {0, 0, 0, lc};
Point(2) = {10, 0, 0, lc};
Point(3) = {10, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

// Define lines
Line(1) = {1, 2}; // Bottom boundary
Line(2) = {2, 3}; // Outlet
Line(3) = {3, 4}; // Top boundary
Line(4) = {4, 1}; // Inlet

// Define line loop and plane surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Define physical groups
surafce_marker = 6;
Physical Line(0) = {4};
Physical Line(1) = {2};
Physical Line(2) = {1, 3};
Physical Surface("surface") = {1};

Mesh 2;
Save "rect.msh";
