//+
rad = DefineNumber[ 1, Name "Parameters/rad" ];
//+
ndens = DefineNumber[ 0.05, Name "Parameters/ndens" ];
//+
Point(1) = {0, 0, 0, ndens};
//+
Point(2) = {rad, 0, 0, ndens};
//+
Point(3) = {-rad, 0, 0, ndens};
//+
Point(4) = {0, rad, 0, ndens};
//+
Point(5) = {0, -rad, 0, ndens};
//+
Point(6) = {0, 0, rad, ndens};
//+
Point(7) = {0, 0, -rad, ndens};
//+
Circle(1) = {2, 1, 4};
//+
Circle(2) = {4, 1, 3};
//+
Circle(3) = {3, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Circle(5) = {2, 1, 6};
//+
Circle(6) = {6, 1, 3};
//+
Circle(7) = {3, 1, 7};
//+
Circle(8) = {7, 1, 2};
//+
Circle(9) = {4, 1, 6};
//+
Circle(10) = {6, 1, 5};
//+
Circle(11) = {5, 1, 7};
//+
Circle(12) = {7, 1, 4};
//+
Line Loop(1) = {2, 7, 12};
//+
Surface(1) = {1};
//+
Line Loop(2) = {2, -6, -9};
//+
Surface(2) = {2};
//+
Line Loop(3) = {3, -10, 6};
//+
Surface(3) = {3};
//+
Line Loop(4) = {3, 11, -7};
//+
Surface(4) = {4};
//+
Line Loop(5) = {4, -8, -11};
//+
Surface(5) = {5};
//+
Line Loop(6) = {4, 5, 10};
//+
Surface(6) = {6};
//+
Line Loop(7) = {1, 9, -5};
//+
Surface(7) = {7};
//+
Line Loop(8) = {1, -12, 8};
//+
Surface(8) = {8};
//+
Surface Loop(1) = {2, 1, 4, 3, 6, 5, 8, 7};
//+
Volume(1) = {1};



// Refine membrane and coarsen cytosol

// Define Ball field
Field[1] = Ball;
Field[1].Radius = 0.95;
//Field[1].Thickness = 0.2; //Only for later versions of gmsh
Field[1].VIn = 0.2;
Field[1].VOut = 0.10;

// Use the Ball field for generating the mesh
Background Field = 1;

// Do not use the boundary cell size in the interior
Mesh.CharacteristicLengthExtendFromBoundary = 0;
