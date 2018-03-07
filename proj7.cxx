/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include "tricase.cxx"
#include "TriangleList.h"

#define ISOVALUE 3.2

// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
    // 2D
    // return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    // return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    idx[0] = pointId%dims[0];
    idx[1] = (pointId/dims[0])%dims[1];
    idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    // idx[0] = pointId%dims[0];
    // idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = cellId/(dims[0]-1);
}

// ****************************************************************************
//  Function: Interp
//
//  Arguments:
//     pt: a one-dimentional location
//     a: location of first point
//     b: location of second point
//     fa: value at first point
//     fb: value at second point
//
//  Returns: the interpolated field value.
// ****************************************************************************

float
Interp(const float pt, const float a, const float b, const float fa, const float fb){
    return (a + (((pt-fa)/(fb-fa))*(b-a)));
}


int
EvaluateIso(int pointIndex, const float *F, const float isovalue) {
    return (F[pointIndex] > isovalue) ? 1:0;
}

void
GetCellPoints(int *points, int cellIndex, const int *dims) {
    
    int idx[3];
    int tmp[3];
    GetLogicalCellIndex(idx, cellIndex, dims);
    
    
    tmp[0] = idx[0];
    tmp[1] = idx[1];
    tmp[2] = idx[2];

    int x, y, z;
    int count = 0;

    for (z = 0; z < 2; z++) {
        tmp[2] += z;
        for (y = 0; y < 2; y++) {
            tmp[1] += y;
            for (x = 0; x < 2; x++) {
                tmp[0] += x;
                points[count] = GetPointIndex(tmp, dims);
                tmp[0] = idx[0];
                count++;
            }
            tmp[1] = idx[1];
        }
    }
}

int
IdentifyCase(int cellIndex, const int *dims, const float *F){

    int point[8];
    GetCellPoints(point, cellIndex, dims);

    int i;
    int total = 0;
    for (i = 0; i < 8; i++) {
        if (EvaluateIso(point[i], F, ISOVALUE))
            total += pow(2, i);
    }

    return total;
}

void
InterpEdge(float *pt, int edgeNum, int cellIndex, const float *X, const float *Y, const float *Z, const float *F, const int *dims) {
    
    int point[8];
    GetCellPoints(point, cellIndex, dims);
    
    static int edges[12][2]=
        { {0,1}, {1,3}, {2,3}, {0,2},
          {4,5}, {5,7}, {6,7}, {4,6},
          {0,4}, {1,5}, {2,6}, {3,7} };

    int idx[2][3];
    int edgePoints[2];
    edgePoints[0] = point[edges[edgeNum][0]];
    edgePoints[1] = point[edges[edgeNum][1]];
    GetLogicalPointIndex(idx[0], edgePoints[0], dims);
    GetLogicalPointIndex(idx[1], edgePoints[1], dims);
    
    pt[0] = Interp(ISOVALUE, X[idx[0][0]], X[idx[1][0]], F[edgePoints[0]], F[edgePoints[1]]);
    pt[1] = Interp(ISOVALUE, Y[idx[0][1]], Y[idx[1][1]], F[edgePoints[0]], F[edgePoints[1]]);
    pt[2] = Interp(ISOVALUE, Z[idx[0][2]], Z[idx[1][2]], F[edgePoints[0]], F[edgePoints[1]]);
    
}

void
EvaluateCells(const int *dims, const float *X, const float *Y, const float *Z, const float *F, TriangleList *tl) {
    int i, j, icase;
    float pt1[3];
    float pt2[3];
    float pt3[3];
        
    for (i = 0; i < GetNumberOfCells(dims); i++) {
        icase = IdentifyCase(i, dims, F);
        int *edges  = triCase[icase];
       
        while (*edges != -1) {

            int v0 = *edges++;
            int v1 = *edges++;
            int v2 = *edges++;
            
            InterpEdge(pt1, v0, i, X, Y, Z, F, dims);
            InterpEdge(pt2, v1, i, X, Y, Z, F, dims);
            InterpEdge(pt3, v2, i, X, Y, Z, F, dims);

            tl->AddTriangle(pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], pt3[0], pt3[1], pt3[2]);
        }
    }
}


int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj6B.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    TriangleList tl;

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!
    EvaluateCells(dims, X, Y, Z, F, &tl);

    vtkPolyData *pd = tl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
