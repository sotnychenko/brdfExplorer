
#include <QtGui>
#include <QCheckBox>
#include "ImageSliceWindow.h"
#include "ParameterWindow.h"
#include "PlotPCASlice2DWindow.h"
#include "BRDFMeasuredMERL.h"

#include <iostream>

using namespace std;

PlotPCASlice2DWindow::PlotPCASlice2DWindow( ParameterWindow* paramWindow )
{
    brdfs = paramWindow->getBRDFList();
    attrListReady=false;

    customPlot = new QCustomPlot(this);
    customPlot->setObjectName(QString::fromUtf8("customPlot"));

    // configure axis rect:
    customPlot->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
    customPlot->axisRect()->setupFullAxesBox(true);
    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("y");



    connect( paramWindow, SIGNAL(brdfListChanged(std::vector<brdfPackage>)), this, SLOT(brdfListChanged(std::vector<brdfPackage>)) );



   customPlot->setMinimumHeight((int)((float)(paramWindow->size().height()/2)*1.55f));
   customPlot->setMinimumWidth(paramWindow->size().width()+(int)((float)(paramWindow->size().width())*.55f));
   customPlot->setMaximumHeight((int)((float)(paramWindow->size().height()/2)*1.55f));
   customPlot->setMaximumWidth(paramWindow->size().width()+(int)((float)(paramWindow->size().width())*.55f));



    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(customPlot);
    mainLayout->setAlignment(Qt::AlignCenter);

    QHBoxLayout *buttonLayout = new QHBoxLayout;
    mainLayout->addLayout(buttonLayout);


    attrComboBox = new QComboBox();

    connect( attrComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(attrComboBoxIndChanged(int)) );

    buttonLayout->addWidget(attrComboBox);

    xAxisComboBox = new QComboBox();
    connect( xAxisComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(xAxisComboBoxIndChanged(int)) );

    buttonLayout->addWidget(xAxisComboBox);

    yAxisComboBox = new QComboBox();
    connect( yAxisComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(yAxisComboBoxIndChanged(int)) );

    buttonLayout->addWidget(yAxisComboBox);

    buttonClean = new QPushButton("Clean", this);

    // add the solo button
    buttonClean->setToolTip( "Cleans the path" );
    connect( buttonClean, SIGNAL(clicked()), this, SLOT(buttonCleanPushed()) );

    buttonClean->setFixedWidth( 50 );
    buttonClean->setFixedHeight( 25 );
    mainLayout->addWidget( buttonClean,0,Qt::AlignRight );

   // mainLayout->setAlignment(Qt::AlignCenter);
    setLayout(mainLayout);

}


void PlotPCASlice2DWindow::buttonCleanPushed()

{
if(brdfs.size()>0){

   brdfs[0].brdf->brdfParam->paths.at(attrComboBox->currentIndex()).alpha.clear();
   drawPlot();
}

return;
}

void PlotPCASlice2DWindow::attrComboBoxIndChanged(int index)

{
    //cout<<"index changed"<<index<<endl;
    attrComboBox->setCurrentIndex( index );
    if(brdfs.size()>0)
    drawPlot();
    return;
}

void PlotPCASlice2DWindow::xAxisComboBoxIndChanged(int index)

{

     xAxisComboBox->setCurrentIndex( index );
     if(brdfs.size()>0)
     drawPlot();
     return;
}
void PlotPCASlice2DWindow::yAxisComboBoxIndChanged(int index)

{
    yAxisComboBox->setCurrentIndex( index );
    if(brdfs.size()>0)
    drawPlot();
    return;
}

void PlotPCASlice2DWindow::brdfListChanged( std::vector<brdfPackage> brdfList )
{
if(brdfList.size()==0)

{
    customPlot->clearPlottables();
  //  customPlot->plotLayout()->remove(colorScale);

    customPlot->replot();

    return;

}

if (!(brdfList[0].brdf->wasProjected)) return;
  // replace the BRDF pointers
  brdfs = brdfList;

  if(!attrListReady)
  {
      for(int i=0;i<brdfs[0].brdf->brdfParam->attrNames.size(); i++)
      attrComboBox->addItem( QString::fromStdString(brdfs[0].brdf->brdfParam->attrNames.at(i)) );

      xAxisComboBox->addItem("alpha1");
      xAxisComboBox->addItem("alpha2");
      xAxisComboBox->addItem("alpha3");
      xAxisComboBox->addItem("alpha4");
      xAxisComboBox->addItem("alpha5");

      yAxisComboBox->addItem("alpha1");
      yAxisComboBox->addItem("alpha2");
      yAxisComboBox->addItem("alpha3");
      yAxisComboBox->addItem("alpha4");
      yAxisComboBox->addItem("alpha5");

      attrComboBox->setCurrentIndex( 0 );
      xAxisComboBox->setCurrentIndex( 0 );
      yAxisComboBox->setCurrentIndex( 1 );
      attrListReady = true;
  }

  drawPlot();


}


void PlotPCASlice2DWindow::drawPlot()
{

    customPlot->clearPlottables();
    customPlot->plotLayout()->remove(colorScale);

//cout<<"path size"<<brdfs[0].brdf->brdfParam->path.x.size()<<endl;
    // set up the QCPColorMap:

    colorMap= new QCPColorMap(customPlot->xAxis, customPlot->yAxis);
    customPlot->addPlottable(colorMap);
    // add a color scale:

    colorScale= new QCPColorScale(customPlot);
    customPlot->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
    colorScale->setType(QCPAxis::atRight);



   // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)

    int nx = 100;
    int ny = 100;

    colorMap->data()->setSize(nx, ny); // we want the color map to have nx * ny data points
    colorMap->data()->setRange(QCPRange(-1.5, 0.5), QCPRange(-1.5, 0.5)); // and span the coordinate range -4..4 in both key (x) and value (y) dimensions
    // now we assign some data, by accessing the QCPColorMapData instance of the color map:
    double x, y, z;


    MatrixXf proj = brdfs[0].brdf->brdfParam->proj;


    MatrixXf proj_Lab =  BRDFMeasuredMERL::rgb2Lab(proj);
    float* xnew = new float[proj_Lab.rows()];

    for (int i = 0; i < proj_Lab.rows(); i++)
        xnew[i] = proj_Lab(i, 0) / 100.0;

    cnpy::NpyArray arr_mv1 = brdfs[0].brdf->brdfParam->npzFiles.at(attrComboBox->currentIndex())["Centers"];
    double* mv1 = reinterpret_cast<double*>(arr_mv1.data);
    Matrix<float, Dynamic, Dynamic, RowMajor> Centers;
    Centers.resize(arr_mv1.shape[1], arr_mv1.shape[0]);
    for (int i = 0; i < arr_mv1.shape[0] * arr_mv1.shape[1]; i++)
        Centers.data()[i] = mv1[i];

    arr_mv1 =  brdfs[0].brdf->brdfParam->npzFiles.at(attrComboBox->currentIndex())["betas"];
    double* betasTemp = reinterpret_cast<double*>(arr_mv1.data);
    float* betas = new float[arr_mv1.shape[0] * arr_mv1.shape[1]];
    for (int i = 0; i < arr_mv1.shape[0] * arr_mv1.shape[1]; i++)
        betas[i] = betasTemp[i];

    arr_mv1 =  brdfs[0].brdf->brdfParam->npzFiles.at(attrComboBox->currentIndex())["Theta"];

    double* ThetaTemp = reinterpret_cast<double*>(arr_mv1.data);
    float* Theta = new float[arr_mv1.shape[0] * arr_mv1.shape[1]];
    for (int i = 0; i < arr_mv1.shape[0] * arr_mv1.shape[1]; i++)
        Theta[i] = ThetaTemp[i];

  //  float ynew = BRDFMeasuredMERL::evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew);

  // cout<<"ynew "<<ynew<<endl;


    for (int xIndex=0; xIndex<nx; ++xIndex)
    {
      for (int yIndex=0; yIndex<ny; ++yIndex)
      {
        colorMap->data()->cellToCoord(xIndex, yIndex, &x, &y);


                      xnew[0]=-0.1582;
                      xnew[1]=-0.0787;
                      xnew[2]=-0.1653;
                      xnew[3]=-0.1170;
                      xnew[4]=-0.0884;

                      xnew[xAxisComboBox->currentIndex()]=x;
                      xnew[yAxisComboBox->currentIndex()]=y;


          z = BRDFMeasuredMERL::evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew);

       if(z>1) z=1.0;
       if(z<0) z =0.0;

        colorMap->data()->setCell(xIndex, yIndex, z);
      }
    }
    colorScale->setDataRange(QCPRange(0.0, 1.0));
    colorMap->setColorScale(colorScale); // associate the color map with the color scale
   // colorScale->axis()->setLabel("Magnetic Field Strength");

    // set the color gradient of the color map to one of the presets:
    colorMap->setGradient(QCPColorGradient::gpThermal);
    // we could have also created a QCPColorGradient instance and added own colors to
    // the gradient, see the documentation of QCPColorGradient for what's possible.

    // rescale the data dimension (color) such that all data points lie in the span visualized by the color gradient:

    //colorMap->rescaleDataRange();


    // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
    QCPMarginGroup *marginGroup = new QCPMarginGroup(customPlot);
    customPlot->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
    colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

    // rescale the key (x) and value (y) axes so the whole color map is visible:

 customPlot->rescaleAxes();
    //customPlot->legend->setVisible(true);
   // customPlot->legend->setFont(QFont("Helvetica", 9));

      // generate data:
    if(brdfs[0].brdf->brdfParam->paths.at(attrComboBox->currentIndex()).alpha.size()!=0)
    {
        QPen pen;
        pen.setWidth(2.5);
          customPlot->addGraph();
          pen.setColor(QColor(0.0, 0.0, 255.0));
          customPlot->graph()->setPen(pen);
          customPlot->graph()->setLineStyle((QCPGraph::LineStyle)1);
          customPlot->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 8.5));


      customPlot->graph()->setData(brdfs[0].brdf->brdfParam->paths.at(attrComboBox->currentIndex()).alpha.at(xAxisComboBox->currentIndex()),
                                   brdfs[0].brdf->brdfParam->paths.at(attrComboBox->currentIndex()).alpha.at(yAxisComboBox->currentIndex()));
      customPlot->graph()->rescaleAxes(true);



    // add the arrow:
    QCPItemLine *arrow = new QCPItemLine(customPlot);
    customPlot->addItem(arrow);


    arrow->start->setCoords(brdfs[0].brdf->brdfParam->paths.at(attrComboBox->currentIndex()).alpha.at(xAxisComboBox->currentIndex()).at(0),
                            brdfs[0].brdf->brdfParam->paths.at(attrComboBox->currentIndex()).alpha.at(yAxisComboBox->currentIndex()).at(0));
    arrow->end->setCoords(brdfs[0].brdf->brdfParam->paths.at(attrComboBox->currentIndex()).alpha.at(xAxisComboBox->currentIndex()).last(),
                          brdfs[0].brdf->brdfParam->paths.at(attrComboBox->currentIndex()).alpha.at(yAxisComboBox->currentIndex()).last()); // point to (4, 1.6) in x-y-plot coordinates
    arrow->setHead(QCPLineEnding::esSpikeArrow);





    }

    customPlot->replot();


}





