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
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <dirent.h>
#include "ParameterGroupWidget.h"
#include "FloatVarWidget.h"
#include "ColorVarWidget.h"
#include "BRDFAnalytic.h"
#include "ParameterWindow.h"
#include "Paths.h"
#include "BRDFMeasuredMERL.h"
#include "projectToPCA.h"
#include"cnpy.h"

using namespace std;

#define NUM_DRAW_COLORS 7
float drawColors[NUM_DRAW_COLORS][3] =
{
    {0.65,0,0},
    {0,0.65,0},
    {0,0,0.65},
    {0.65,0.65,0},
    {0,0.65,0.65},
    {0.65,0,0.65},
    {0.65,0.65,0.65}
};

std::string extractFilename( const std::string& path )
{
    return path.substr( path.find_last_of( '/' ) +1 );
}

int ParameterGroupWidget::colorIndex = 0;

ParameterGroupWidget::ParameterGroupWidget( ParameterWindow* pWindow, BRDFBase* brdfIn )
    : QFrame(), paramWindow(pWindow), brdf(brdfIn), visibleCheckBox(NULL), dirty(true)
{
    // start by setting the draw color of the BRDF
    brdfDrawColor[0] = drawColors[colorIndex][0];
    brdfDrawColor[1] = drawColors[colorIndex][1];
    brdfDrawColor[2] = drawColors[colorIndex][2];

    colorIndex = (colorIndex + 1) % NUM_DRAW_COLORS;

    dosomething =false;

    // now let's get to the layout
    QVBoxLayout *layout = new QVBoxLayout;
    layout->setMargin( 0 );
    layout->setContentsMargins( 0, 0, 0, 0 );

    // the parameter window needs to know to emit changes whenever this BRDF is reloaded
    connect( this, SIGNAL(brdfChanged(ParameterGroupWidget*)), paramWindow, SLOT(emitBRDFListChanged()) );


    // add and connect the show/hide button
    titleButton = new QPushButton( QString("  ") + QString(extractFilename(brdf->getName()).c_str()) );
    titleButton->setFixedHeight( 22 ); 
    connect( titleButton, SIGNAL(clicked()), this, SLOT(titleButtonPushed()) );


    // set the button color to the BRDF's draw color
    changeTitleButtonColor( false );



    // make the command QFrame - contains the visible checkbox and reload/close buttons
    QFrame* cmdFrame = new QFrame;

    QHBoxLayout* cmdLayout = new QHBoxLayout;
    cmdLayout->setMargin( 0 );
    cmdLayout->setContentsMargins( 0, 0, 0, 0 );
    cmdLayout->setSpacing( 2 );
    cmdFrame->setLayout( cmdLayout );


    visibleCheckBox = new QCheckBox( "Visible" );
    visibleCheckBox->setChecked( true );
    cmdLayout->addWidget( visibleCheckBox );
    connect( visibleCheckBox, SIGNAL(stateChanged(int)), this, SLOT(paramChanged()) );
    
   std::string extension =  brdfIn->getName().substr( brdfIn->getName().find_last_of( '.' ) +1 );
    pcaCheckBox=NULL;
    if( extension == "binary" ){
    pcaCheckBox = new QCheckBox( "Project to PCA" );
    pcaCheckBox->setChecked( false );
    cmdLayout->addWidget( pcaCheckBox );
    connect( pcaCheckBox, SIGNAL(stateChanged(int)), this, SLOT(paramChanged()) );
    
   // colorSpaceBox = new QCheckBox( "2nd ver" );
   // colorSpaceBox->setChecked( false );
  //  cmdLayout->addWidget( colorSpaceBox );
   // connect( colorSpaceBox, SIGNAL(stateChanged(int)), this, SLOT(paramChanged()) );
    buttonReset = new QPushButton("reset", this);
    // add the solo button
    buttonReset->setToolTip( "resets the brdf projection" );
    connect( buttonReset, SIGNAL(clicked()), this, SLOT(buttonResetPushed()) );

    buttonReset->setFixedWidth( 35 );
    buttonReset->setFixedHeight( 25 );
    cmdLayout->addWidget( buttonReset,0,Qt::AlignRight );

   }


    // add the solo button
    soloButton = new QPushButton();
    QPixmap* soloPixmap = new QPixmap((getImagesPath() + "soloSmall.png").c_str());
    soloButton->setIconSize( QSize(soloPixmap->width(), soloPixmap->height()) );
    soloButton->setIcon( QIcon(*soloPixmap) );
    soloButton->setFixedWidth( 30 );
    soloButton->setFixedHeight( 24 );
    soloButton->setCheckable( true );
    soloButton->setChecked( false );
    soloButton->setToolTip( "Solo this BRDF" );
    connect( soloButton, SIGNAL(clicked()), this, SLOT(soloButtonPushed()) );
    cmdLayout->addWidget( soloButton );


    // add the solo colors button
    soloColorsButton = new QPushButton();
    QPixmap* soloColorsPixmap = new QPixmap((getImagesPath() + "soloColorsSmall.png").c_str());
    soloColorsButton->setIconSize( QSize(soloColorsPixmap->width(), soloColorsPixmap->height()) );
    soloColorsButton->setIcon( QIcon(*soloColorsPixmap) );
    soloColorsButton->setFixedWidth( 30 );
    soloColorsButton->setFixedHeight( 24 );
    soloColorsButton->setCheckable( true );
    soloColorsButton->setChecked( false );
    soloColorsButton->setToolTip( "Solo this BRDF's color channels" );
    connect( soloColorsButton, SIGNAL(clicked()), this, SLOT(soloColorsButtonPushed()) );
    cmdLayout->addWidget( soloColorsButton );
   
    
    QMenu* optionsMenu = new QMenu;
    QAction* reloadAction = optionsMenu->addAction( "Reload BRDF" );
    QPixmap* reloadPixmap = new QPixmap((getImagesPath() + "reloadSmall.png").c_str());
    reloadAction->setIcon( QIcon(*reloadPixmap) );    
    connect( reloadAction, SIGNAL(triggered()), this, SLOT(reloadButtonPushed()) );
    
    QAction* resetAction = optionsMenu->addAction( "Reload BRDF and reset to default" );
    QPixmap* resetPixmap = new QPixmap((getImagesPath() + "resetSmall.png").c_str());
    resetAction->setIcon( QIcon(*resetPixmap) );
    connect( resetAction, SIGNAL(triggered()), this, SLOT(resetButtonPushed()) );


    QAction* saveBRDFAction = optionsMenu->addAction( "Save BRDF..." );
    QPixmap* folderPixmap = new QPixmap((getImagesPath() + "folderSmall.png").c_str());
    saveBRDFAction->setIcon( QIcon(*folderPixmap) );
    connect( saveBRDFAction, SIGNAL(triggered()), this, SLOT(saveBRDFButtonPushed()) );
    
    QAction* saveAction = optionsMenu->addAction( "Save Parameters File..." );

    saveAction->setIcon( QIcon(*folderPixmap) );
    connect( saveAction, SIGNAL(triggered()), this, SLOT(saveParamsFileButtonPushed()) );
    
    optionsMenu->addSeparator();
    
    QAction* closeAction = optionsMenu->addAction( "Close BRDF" );
    QPixmap* closePixmap = new QPixmap((getImagesPath() + "closeSmall.png").c_str());
    closeAction->setIcon( QIcon(*closePixmap) );
    connect( closeAction, SIGNAL(triggered()), this, SLOT(removeButtonPushed()) );
    

    // add the button with the menu dropdown
    QPushButton* menuButton = new QPushButton();
    menuButton->setFixedWidth( 24 );
    menuButton->setFixedHeight( 24 );
    menuButton->setMenu( optionsMenu );    
    cmdLayout->addWidget( menuButton );


    // make the container frame and its layout
    containerFrame = new QFrame;
    containerLayout = new QVBoxLayout( 0 );
    containerLayout->setMargin( 0 );
    containerLayout->setContentsMargins( 0, 0, 0, 0 );
    containerLayout->addWidget( cmdFrame );
    containerLayout->setSpacing( 2 );
    containerFrame->setLayout( containerLayout );

    addParameterWidgets();

    // add the widgets to the layout
    layout->addWidget( titleButton );
    layout->addWidget( containerFrame );


    setLayout(layout);
}


ParameterGroupWidget::~ParameterGroupWidget()
{
    delete brdf;
}

void ParameterGroupWidget::buttonResetPushed()

{




     BRDFMeasuredMERL* mb = new BRDFMeasuredMERL;

     mb->brdfParam=brdf->brdfParam;
  mb->reset( brdf->numBRDFSamples,brdf->getName().c_str());
       brdf->wasProjected =false;
   updateAttrSliders( mb->brdfParam);

delete brdf;
brdf=mb;
brdf->numBRDFSamples=mb->numBRDFSamples;
brdf->brdfParam=mb->brdfParam;
brdf->wasProjected =true;

for(int i=0; i<brdf->brdfParam->paths.size();i++)
 brdf->brdfParam->paths.at(i).alpha.clear();

paramChanged();

}

void ParameterGroupWidget::titleButtonPushed()
{
	// toggle the visibility of the parameters in the container frame
	bool visible = containerFrame->isVisible();
	containerFrame->setVisible( !visible );
}


void ParameterGroupWidget::removeButtonPushed()
{
	emit( removeBRDFButtonPushed(this) );
}
void ParameterGroupWidget::addAttributeWidgets()
{
 FloatVarWidget* fv;
    
      DIR *dir;
      struct dirent *ent;
      int it=0;
      if ((dir = opendir ("trained_RBFN_npz\\")) != NULL) {
         /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
       
    
          std::string str(ent->d_name);
    
           if(str.substr(str.find_last_of(".") + 1) == "npz") {
      
   
    size_t found = str.find("Name_");
    
    if(found!=std::string::npos){ 
     
     string strNew = str.substr (found,str.size()-1);
   
     unsigned first = strNew.find("_");
     strNew=  strNew.substr (first+1,strNew.size()-1);
     unsigned last = strNew.find("_");
     brdf->brdfParam->attrNames.push_back(strNew.substr (0,last));
       
       fv = new FloatVarWidget(QString::fromStdString(strNew.substr (0,last)), 0.2, 1.0, 1.0);
    }
    else    
    fv = new FloatVarWidget(QString::number(it), 0.2, 1.0, 1.0);
    fv->setId(it);


    
    connect(fv, SIGNAL(valueChanged(float,int)), this, SLOT(attrChanged(float,int)));
    
    this->layout()->addWidget(fv);
  
    cnpy::npz_t my_npz = cnpy::npz_load("trained_RBFN_npz\\"+str);
   brdf->brdfParam->npzFiles.push_back(my_npz); 

 //   printf ("%s\n", ent->d_name);
    it++;
  }
  
  }
  closedir (dir);
        } else {
  /* could not open directory */
  perror ("");

       }
       
 
  }
  void  ParameterGroupWidget::updateAttrSliders(brdfMERLparam* brdfParam)
  {
 
  int offset = this->layout()->count()-brdfParam->npzFiles.size();
        for(int i=0;i<brdfParam->npzFiles.size();i++)
      {  
            FloatVarWidget* w = (FloatVarWidget*)(this->layout()->itemAt(i+offset)->widget());      
             w->setValue( brdfParam->attrValues[i]);


       }
  }

BRDFBase* ParameterGroupWidget::getUpdatedBRDF()
{
 // if this isn't visible, then don't send back a BRDF
    if( visibleCheckBox && !visibleCheckBox->isChecked() && !isSoloing() )
        return NULL;

if(dosomething==true)
{


        dosomething=false;

         BRDFMeasuredMERL* mb = new BRDFMeasuredMERL;

         mb->brdfParam=brdf->brdfParam;
      mb->projectShort( brdf->numBRDFSamples,brdf->getName().c_str());
           brdf->wasProjected =false;
       updateAttrSliders( mb->brdfParam);

  delete brdf;
 brdf=mb;
 brdf->numBRDFSamples=mb->numBRDFSamples;
brdf->brdfParam=mb->brdfParam;
 brdf->wasProjected =true;


}

if(pcaCheckBox!=NULL)
if( pcaCheckBox->isChecked() && !brdf->wasProjected )

  {

    std::cerr<<brdf->getName().c_str()<<std::endl;
    brdf->brdfParam =new brdfMERLparam;
   // brdf->brdfParam->verOfColorSpace =  colorSpaceBox->isChecked();
    addAttributeWidgets();

      BRDFMeasuredMERL* mb = new BRDFMeasuredMERL;
      mb->brdfParam= brdf->brdfParam;
    mb->project(brdf->getName().c_str());

    updateAttrSliders( mb->brdfParam);

delete brdf;
brdf=mb;
brdf->numBRDFSamples=mb->numBRDFSamples;
brdf->brdfParam=mb->brdfParam;

      brdf->wasProjected =true;

  }


    int floatIndex = 0;
    int boolIndex = 0;
	int colorIndex = 0;

    for( int i = 0; i < (int)brdfParamWidgets.size(); i++ )
    {
        // float
        if( brdfParamWidgets[i].type == BRDF_VAR_FLOAT )
        {
            float value = 0.0;
            FloatVarWidget* w = dynamic_cast<FloatVarWidget*>(brdfParamWidgets[i].widget);
            if( w ) value = w->getValue();
            brdf->setFloatParameterValue( floatIndex++, value );
        }

        // bool
        else if( brdfParamWidgets[i].type == BRDF_VAR_BOOL )
        {
            float value = 0.0;
            QCheckBox* c = dynamic_cast<QCheckBox*>(brdfParamWidgets[i].widget);
            if( c ) value = c->isChecked() ? 1.0 : 0.0;
            brdf->setBoolParameterValue( boolIndex++, value );
            
        }

        // colors
        else if( brdfParamWidgets[i].type == BRDF_VAR_COLOR )
        {
			QColor value = QColor(1,1,1);
            ColorVarWidget* c = dynamic_cast<ColorVarWidget*>(brdfParamWidgets[i].widget);
            if( c ) value = c->getValue();
            brdf->setColorParameterValue( colorIndex++, value.redF(), value.greenF(), value.blueF() );  
        }
    }

    return brdf;
}


void ParameterGroupWidget::reloadButtonPushed()
{
    reload(false);
}

void ParameterGroupWidget::resetButtonPushed()
{
    reload(true);
}
void ParameterGroupWidget::attrChanged(float v,int id)
{
 if(! brdf->wasProjected) return;
dosomething=true;
brdf->brdfParam->newAttrVal = v;
brdf->brdfParam->idOfVal=id;

paramChanged();

}

void ParameterGroupWidget::reload(bool resetToDefaults)
{
    // attempt to create a new BRDF based on this one
    BRDFBase* newBRDF = brdf->cloneBRDF(resetToDefaults);
    if( !newBRDF )
            return;

    // remove all the items from the layout
    QLayoutItem *child;
    while ((child = containerLayout->takeAt(1)) != 0)
    {
        child->widget()->deleteLater();
        delete child;
    }
    brdfParamWidgets.clear();

    // switch out the old BRDF for the new one
    delete brdf;
    brdf = newBRDF;

    // now add parameter widgets based on the new BRDF
    addParameterWidgets();

    // now redraw
    paramChanged();
}


void ParameterGroupWidget::paramChanged()
{
    setDirty( true );
    emit( brdfChanged(this) );
}


void ParameterGroupWidget::addParameterWidgets()
{
    // no BRDF? No parameters to load.
    if( !brdf )
            return;

    // loop through and add all the parameters to the list
    int float_index = 0, bool_index = 0, color_index = 0;
    for( int i = 0; i < brdf->getParameterCount(); i++ )
    {
        switch (brdf->getParameterType(i)) {
            case BRDF_VAR_FLOAT:
            {
                const brdfFloatParam* p = brdf->getFloatParameter( float_index++ );
                FloatVarWidget* fv = new FloatVarWidget(QString(p->name.c_str()), p->minVal, p->maxVal, p->defaultVal);
                fv->setValue( p->currentVal );
                connect(fv, SIGNAL(valueChanged(float)), this, SLOT(paramChanged()));
                brdfParamWidgets.push_back( BRDFParamWidget(BRDF_VAR_FLOAT, fv) );
                containerLayout->addWidget( fv );
            }
            break;

            case BRDF_VAR_BOOL:
            {
                const brdfBoolParam* p = brdf->getBoolParameter( bool_index++ );
                QCheckBox* bv = new QCheckBox(QString(p->name.c_str()));
                bv->setChecked( bool(p->currentVal > 0.0001) );
                connect(bv, SIGNAL(stateChanged(int)), this, SLOT(paramChanged()));
                brdfParamWidgets.push_back( BRDFParamWidget(BRDF_VAR_BOOL, bv) );
                containerLayout->addWidget( bv );
            }
            break;

            case BRDF_VAR_COLOR:
            {
                const brdfColorParam* p = brdf->getColorParameter( color_index++ );
                ColorVarWidget* cv = new ColorVarWidget(QString(p->name.c_str()),
                                                        p->currentVal[0], p->currentVal[1], p->currentVal[2]);
                connect(cv, SIGNAL(valueChanged()), this, SLOT(paramChanged()));
                brdfParamWidgets.push_back( BRDFParamWidget(BRDF_VAR_COLOR, cv) );
                containerLayout->addWidget( cv );
            }
            break;
        }
    }
}


void ParameterGroupWidget::brdfBeingSoloed( ParameterGroupWidget* pgw, bool withColors )
{
    // is this the BRDF being soloed?
    if( pgw == this )
    {
        if( withColors )
        {
            soloColorsButton->setChecked( true );
            soloButton->setChecked( false );
        }
        else
        {
            soloButton->setChecked( true );
            soloColorsButton->setChecked( false );
        }
        changeTitleButtonColor( false );
    }

    // no BRDF being soloed
    else if( pgw == NULL )
    {
        changeTitleButtonColor( false );
        soloColorsButton->setChecked( false );
        soloButton->setChecked( false );
    }

    // some other BRDF being soloed
    else
    {
        soloColorsButton->setChecked( false );
        soloButton->setChecked( false );
        changeTitleButtonColor( true );
    }

}



bool ParameterGroupWidget::isSoloing()
{
    if( !soloColorsButton )
        return false;
    if( !soloButton )
        return false;

    return soloColorsButton->isChecked() || soloButton->isChecked();
}



void ParameterGroupWidget::soloButtonPushed()
{
    // is solo mode for this BRDF being turned on?
    if( isSoloing() )
    {
        //visibleCheckBox->setChecked( true );
        emit( enteringSoloMode(this, false) );
    }

    // nope, we're turning it on
    else
    {
        emit( enteringSoloMode(NULL, false) );
    }
}


void ParameterGroupWidget::soloColorsButtonPushed()
{
    // is solo mode for this BRDF being turned on?
    if( isSoloing() )
    {
        //visibleCheckBox->setChecked( true );
        emit( enteringSoloMode(this, true) );
    }

    // nope, we're turning it on
    else
    {
        emit( enteringSoloMode(NULL, false) );
    }
}

void ParameterGroupWidget::saveBRDFButtonPushed()
{
    // sanity check - shouldn't happen

    if( !brdf )
        return;
   if(!brdf->wasProjected)
        return;

    // need to come up with a default name for this save file
    QString shortBRDFName = titleButton->text().trimmed();
    int dotIndex = shortBRDFName.lastIndexOf( "." );
    if( dotIndex != -1 )
        shortBRDFName = shortBRDFName.left( dotIndex );


    // get a save file name
    QString fileName = QFileDialog::getSaveFileName(this, "Save BRDF binary File",
                                                          QString("./") + shortBRDFName + QString(".binary"),
                                                          "MERL BRDF (*.binary)" );

    // if we got a filename back... save it
    if( fileName.length() )
        brdf->saveBRDF( fileName.toAscii().constData() );
}
void ParameterGroupWidget::saveParamsFileButtonPushed()
{   
    // sanity check - shouldn't happen
    if( !brdf )
        return;
    
    
    // need to come up with a default name for this save file
    QString shortBRDFName = titleButton->text().trimmed();
    int dotIndex = shortBRDFName.lastIndexOf( "." );
    if( dotIndex != -1 )
        shortBRDFName = shortBRDFName.left( dotIndex );
    

    // get a save file name
    QString fileName = QFileDialog::getSaveFileName(this, "Save BRDF Parameters File",
                                                          QString("./") + shortBRDFName + QString(".bparam"),
                                                          "Parameter Files (*.bparam)" );
    
    // if we got a filename back... save it
    if( fileName.length() )
        brdf->saveParamsFile( fileName.toAscii().constData() );
}


void ParameterGroupWidget::changeTitleButtonColor( bool black )
{
    if( black )
    {
        titleButton->setStyleSheet( "background-color: black; color:white; text-align:left;" );
    }

    else
    {
        char stylesheet[255];
        sprintf( stylesheet, "background-color: #%02x%02x%02x; text-align:left; color:white",
                 int(brdfDrawColor[0]*255.0f), int(brdfDrawColor[1]*255.0f), int(brdfDrawColor[2]*255.0f) );
        titleButton->setStyleSheet( QString(stylesheet) );
    }
}


QColor ParameterGroupWidget::getDrawColor()
{
    QColor c;
    c.setRgbF( brdfDrawColor[0], brdfDrawColor[1], brdfDrawColor[2] );
    return c;
}



