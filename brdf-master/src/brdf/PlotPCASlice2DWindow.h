#ifndef PLOTPCASLICE2DWINDOW_H
#define PLOTPCASLICE2DWINDOW_H



#include <QWidget>
#include "BRDFBase.h"
#include "ShowingBase.h"
#include "qcustomplot.h"

class ParameterWindow;

class QComboBox;


class PlotPCASlice2DWindow : public QWidget
{
    Q_OBJECT


public:
    PlotPCASlice2DWindow( ParameterWindow* pm);

public slots:
    void brdfListChanged( std::vector<brdfPackage> );
    void attrComboBoxIndChanged(int);
    void xAxisComboBoxIndChanged(int);
    void yAxisComboBoxIndChanged(int);

private:
    QCustomPlot* customPlot;
    QComboBox* attrComboBox;
    QComboBox* xAxisComboBox;
    QComboBox* yAxisComboBox;
    QCPColorMap *colorMap;
    QCPColorScale *colorScale;
    std::vector<brdfPackage> brdfs;
    bool attrListReady;
    void drawPlot();
};

#endif // PLOTPCASLICE2DWINDOW_H
