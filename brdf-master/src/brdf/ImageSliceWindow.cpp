/*
Copyright Disney Enterprises, Inc. All rights reserved.

This license governs use of the accompanying software. If you use the software, you
accept this license. If you do not accept the license, do not use the software.

1. Definitions
The terms "reproduce," "reproduction," "derivative works," and "distribution" have
the same meaning here as under U.S. copyright law. A "contribution" is the original
software, or any additions or changes to the software. A "contributor" is any person
that distributes its contribution under this license. "Licensed patents" are a
contributor's patent claims that read directly on its contribution.

2. Grant of Rights
(A) Copyright Grant- Subject to the terms of this license, including the license
conditions and limitations in section 3, each contributor grants you a non-exclusive,
worldwide, royalty-free copyright license to reproduce its contribution, prepare
derivative works of its contribution, and distribute its contribution or any derivative
works that you create.
(B) Patent Grant- Subject to the terms of this license, including the license
conditions and limitations in section 3, each contributor grants you a non-exclusive,
worldwide, royalty-free license under its licensed patents to make, have made,
use, sell, offer for sale, import, and/or otherwise dispose of its contribution in the
software or derivative works of the contribution in the software.

3. Conditions and Limitations
(A) No Trademark License- This license does not grant you rights to use any
contributors' name, logo, or trademarks.
(B) If you bring a patent claim against any contributor over patents that you claim
are infringed by the software, your patent license from such contributor to the
software ends automatically.
(C) If you distribute any portion of the software, you must retain all copyright,
patent, trademark, and attribution notices that are present in the software.
(D) If you distribute any portion of the software in source code form, you may do
so only under this license by including a complete copy of this license with your
distribution. If you distribute any portion of the software in compiled or object code
form, you may only do so under a license that complies with this license.
(E) The software is licensed "as-is." You bear the risk of using it. The contributors
give no express warranties, guarantees or conditions. You may have additional
consumer rights under your local laws which this license cannot change.
To the extent permitted under your local laws, the contributors exclude the
implied warranties of merchantability, fitness for a particular purpose and non-
infringement.
*/

#include <QtGui>
#include <QCheckBox>
#include "ImageSliceWindow.h"
#include "ParameterWindow.h"
#include "FloatVarWidget.h"
#include "ImageSliceWidget.h"

ImageSliceWindow::ImageSliceWindow( ParameterWindow* paramWindow )
{
    glWidget = new ImageSliceWidget( this, paramWindow->getBRDFList() );

    connect( paramWindow, SIGNAL(incidentDirectionChanged(float,float)), glWidget, SLOT(incidentDirectionChanged(float,float)) );
    connect( paramWindow, SIGNAL(brdfListChanged(std::vector<brdfPackage>)), glWidget, SLOT(brdfListChanged(std::vector<brdfPackage>)) );
    
    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(glWidget);

    QHBoxLayout *buttonLayout = new QHBoxLayout;
    mainLayout->addLayout(buttonLayout);

    useThetaHSquared = new QCheckBox( "Square ThetaH" );
    useThetaHSquared->setChecked( false );
    connect( useThetaHSquared, SIGNAL(stateChanged(int)), glWidget, SLOT(useThetaHSquaredChanged(int)) );
    buttonLayout->addWidget(useThetaHSquared);

    showChroma = new QCheckBox( "Show Chroma" );
    showChroma->setChecked( false );
    connect( showChroma, SIGNAL(stateChanged(int)), glWidget, SLOT(showChromaChanged(int)) );
    buttonLayout->addWidget(showChroma);

    FloatVarWidget* fv;
    
    fv = new FloatVarWidget("PhiD", 0, 180.0, 90.0);
    connect(fv, SIGNAL(valueChanged(float)), glWidget, SLOT(phiDChanged(float)));
    mainLayout->addWidget(fv);
    
    fv = new FloatVarWidget("Gamma", 1.0, 5.0, 2.2);
    connect(fv, SIGNAL(valueChanged(float)), glWidget, SLOT(gammaChanged(float)));
    mainLayout->addWidget(fv);
    
    fv = new FloatVarWidget("Exposure", -6.0, 6.0, 0.0);
    connect(fv, SIGNAL(valueChanged(float)), glWidget, SLOT(exposureChanged(float)));
    mainLayout->addWidget(fv);


    setLayout(mainLayout);
    
    setWindowTitle( "Lit Sphere" );
}


void ImageSliceWindow::setShowing( bool s )
{
    if( glWidget )
        glWidget->setShowing( s );
}
