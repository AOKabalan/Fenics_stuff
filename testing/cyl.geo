D = 1.0;
R = 0.5*D;

CX = 4.5*D;
CY = 0.0;

S = 1.0/Sqrt(2);

DR1 = R*S;
DR2 = (R + 2.0*D)*S;

N1 = 35;
N2 = 90;

Point(1) = {0, -4.5*D, 0};
Point(2) = {0, CY - DR2, 0};
Point(3) = (0, CY + DR2, 0);
Point(4) = {0, 4.5*D, 0};
Point(5) = {CX - DR2, -4.5*D, 0};
Point(6) = {CX - DR2, CY - DR2, 0};
Point(7) = {CX - DR2, CY + DR2, 0};
Point(8) = {CX - DR2, 4.5*D, 0};
Point(9) = {CX - DR1, CY - DR1, 0};
Point(10) = {CX - DR1, CY + DR1, 0};
Point(11) = {CX, CY, 0};
Point(12) = {CX + DR1, CY - DR1, 0};
Point(13) = {CX + DR1, CY + DR1, 0};
Point(14) = {CX + DR2, -4.5*D, 0};
Point(15) = {CX + DR2, CY - DR2, 0};
Point(16) = {CX + DR2, CY + DR2, 0};
Point(17) = {CX + DR2, 4.5*D, 0};
Point(18) = {21*D, -4.5*D, 0};
Point(19) = {21*D, CY - DR2, 0};
Point(20) = {21*D, CY + DR2, 0};
Point(21) = {21*D, 4.5*D, 0};
Point(22) = {0, CY, 0};
Point(23) = {CX - 2*D, CY, 0};
Point(24) = {CX + 2*D, CY, 0};
Point(25) = {21*D, CY, 0};
Point(26) = {CX, -4.5*D, 0};
Point(27) = {CX, CY - 2*D, 0};
Point(28) = {CX, CY + 2*D, 0};
Point(29) = {CX, 4.5*D, 0};
Point(30) = {CX, CY - R, 0};
Point(31) = {CX - R, CY, 0};
Point(32) = {CX, CY + R, 0};
Point(33) = {CX + R, CY, 0};
Line(1) = {4, 8};
Line(2) = {8, 29};
Line(3) = {29, 17};
Line(4) = {17, 21};
Line(5) = {21, 20};
Line(6) = {20, 25};
Line(7) = {25, 19};
Line(8) = {19, 18};
Line(9) = {18, 14};
Line(10) = {14, 26};
Line(11) = {26, 5};
Line(12) = {5, 1};
Line(13) = {1, 2};
Line(14) = {2, 22};
Line(15) = {22, 3};
Line(16) = {3, 4};
Line(17) = {3, 7};
Line(18) = {7, 8};
Line(19) = {2, 6};
Line(20) = {6, 5};
Line(21) = {15, 14};
Line(22) = {15, 19};
Line(23) = {16, 17};
Line(24) = {16, 20};
Line(25) = {25, 24};
Line(26) = {23, 22};
Line(27) = {28, 29};
Line(28) = {27, 26};
Line(29) = {27, 15};
Line(30) = {15, 24};
Line(31) = {24, 16};
Line(32) = {16, 28};
Line(33) = {28, 7};
Line(34) = {7, 23};
Line(35) = {23, 6};
Line(36) = {6, 27};
Line(37) = {30, 27};
Line(38) = {31, 23};
Line(39) = {33, 24};
Line(40) = {32, 28};
Line(41) = {10, 7};
Line(42) = {9, 6};
Line(43) = {13, 16};
Line(44) = {12, 15};
Circle(45) = {31, 11, 10};
Circle(46) = {10, 11, 32};
Circle(47) = {32, 11, 13};
Circle(48) = {13, 11, 33};
Circle(49) = {33, 11, 12};
Circle(50) = {12, 11, 30};
Circle(51) = {30, 11, 9};
Circle(52) = {9, 11, 31};
Line Loop(53) = {17, 18, -1, -16};
Ruled Surface(54) = {53};
Line Loop(55) = {33, 18, 2, -27};
Ruled Surface(56) = {55};
Line Loop(57) = {32, 27, 3, -23};
Ruled Surface(58) = {57};
Line Loop(59) = {24, -5, -4, -23};
Ruled Surface(60) = {59};
Line Loop(61) = {25, 31, 24, 6};
Ruled Surface(62) = {61};
Line Loop(63) = {30, -25, 7, -22};
Ruled Surface(64) = {63};
Line Loop(65) = {8, 9, -21, 22};
Ruled Surface(66) = {65};
Line Loop(67) = {21, 10, -28, 29};
Ruled Surface(68) = {67};
Line Loop(69) = {11, -20, 36, 28};
Ruled Surface(70) = {69};
Line Loop(71) = {12, 13, 19, 20};
Ruled Surface(72) = {71};
Line Loop(73) = {19, -35, 26, -14};
Ruled Surface(74) = {73};
Line Loop(75) = {34, 26, 15, 17};
Ruled Surface(76) = {75};
Line Loop(77) = {42, -35, -38, -52};
Ruled Surface(78) = {77};
Line Loop(79) = {36, -37, 51, 42};
Ruled Surface(80) = {79};
Line Loop(81) = {29, -44, 50, 37};
Ruled Surface(82) = {81};
Line Loop(83) = {44, 30, -39, 49};
Ruled Surface(84) = {83};
Line Loop(85) = {39, 31, -43, 48};
Ruled Surface(86) = {85};
Line Loop(87) = {43, 32, -40, 47};
Ruled Surface(88) = {87};
Line Loop(89) = {46, 40, 33, -41};
Ruled Surface(90) = {89};
Line Loop(91) = {38, -34, -41, -45};
Ruled Surface(92) = {91};

RHO1 = 40; // Vertial near-wall parts
RHO2 = 30;
RHO3 = 30;
RHO4 = 200; // Outlet channel
RHO5 = 40;
RHO6 = 40; // Inlet channel

Transfinite Line {13, 20, 28, 21, 8} = RHO1/1.4;
Transfinite Line {16, 18, 27, 23, 5} = RHO1/1.4;
Transfinite Line {15, 34, 45, 48, 31, 6} = RHO2/1.4;
Transfinite Line {14, 35, 52, 49, 30, 7} = RHO2/1.4;
Transfinite Line {11, 36, 51, 46, 33, 2} = RHO3/1.4;
Transfinite Line {10, 29, 50, 47, 32, 3} = RHO3/1.4;
Transfinite Line {9, 22, 25, 24, 4} = RHO4/1.4;
Transfinite Line {42, 37, 44, 39, 43, 40, 41, 38} = RHO5/1.4 Using Progression 1.05;
Transfinite Line {1, 17, 26, 19, 12} = RHO6/1.4;

Transfinite Surface "*";
//Recombine Surface "*"; uncomment for Q elements


//Inlet
Physical Curve(1) = {13,14,15,16}; 
//Outlet
Physical Curve(2) = {5, 6, 7, 8};
//Walls
Physical Curve(3) = {1, 2, 3, 4, 9, 10, 11, 12};
//Cylinder
Physical Curve(4) = {45, 46, 47, 48, 49, 50, 51, 52};
//My surface
Physical Surface(5) = {54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92};
Mesh 2;
Save "flowpastcylinder.msh";
