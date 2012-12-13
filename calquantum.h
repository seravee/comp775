#include <vtkSmartPointer.h>

#include <vtksrep.h>
#include <vtksrepinterpolatemedialsheet.h>
#include <vtksrepvisuprimitives.h>

#include <vtkCellArray.h>
#include <vtkQuad.h>
#include <vtkLine.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>

#include <vtksrepinterpolatemedialcrestcurve.h>
#include <vtksrepinterpolatemedialspokes.h>
#include <vtksrepinterpolatemedialspokeshermite.h>
#include <vtksrepinterpolatecrestspokes.h>
#include <vtksrepinterpolatecrestspokesquartic.h>



void calEigen(double matrix[][2], vector<double> data1, vector<double> data2, int datanum, double &eigen1, double &eigen2);

double calquantum(vtkSmartPointer<vtkSRep> srepfig, vtkSmartPointer<vtkSRepInterpolateMedialSheet> medialsheetinterpolator, vtkSmartPointer<vtkSRep> srepcrest, vtkSmartPointer<vtkSRepInterpolateMedialSpokesHermite> medialspokesinterpolator);
