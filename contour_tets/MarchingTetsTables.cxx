#include "MarchingTetsTables.h"

const int MarchingTetsTables::edgeVertices[6][2] = { {0, 1}, {1, 2}, {2, 0},
                                                     {0, 3}, {1, 3}, {2, 3} };

const int MarchingTetsTables::caseTrianglesEdges[16][7] = {
  {-1, -1, -1, -1, -1, -1, -1},
  { 3,  0,  2, -1, -1, -1, -1},
  { 1,  0,  4, -1, -1, -1, -1},
  { 2,  3,  4,  2,  4,  1, -1},
  { 2,  1,  5, -1, -1, -1, -1},
  { 5,  3,  1,  1,  3,  0, -1},
  { 2,  0,  5,  5,  0,  4, -1},
  { 5,  3,  4, -1, -1, -1, -1},
  { 4,  3,  5, -1, -1, -1, -1},
  { 4,  0,  5,  5,  0,  2, -1},
  { 5,  0,  3,  1,  0,  5, -1},
  { 2,  5,  1, -1, -1, -1, -1},
  { 4,  3,  1,  1,  3,  2, -1},
  { 4,  0,  1, -1, -1, -1, -1},
  { 2,  0,  3, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1}
};

const unsigned MarchingTetsTables::numberOfTriangles[16] = { 0, 1, 1, 2,
                                                             1, 2, 2, 1,
                                                             1, 2, 2, 1,
                                                             2, 1, 1, 0 };

