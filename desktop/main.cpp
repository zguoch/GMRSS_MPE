
// QT includes 
#include <QApplication>
#include <QSurfaceFormat>

#include "QVTKOpenGLWidget.h"
#include "mainwindow.h"

#include <QApplication>

#include "GMRSS_MPEVersion.h"
#include "GMRSS_MPEbash.h"

extern int qInitResources_icons();

int main(int argc, char *argv[])
{
   #ifdef _WIN32
      if(argc == 1)
      {
        // needed to ensure appropriate OpenGL context is created for VTK rendering.
            QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());   // must be here
              // QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
            // QT Stuff
            QApplication app( argc, argv );

            QApplication::setStyle("fusion");

            qInitResources_icons();

            MainWindow myMainWindow;
            //  myMainWindow.showMaximized();
            myMainWindow.show();

            return app.exec();
      }else
      {
        GMRSS_MPEbash::bash_run(argc, argv);
      }
  #else
      if(argc==2 || argc==1)
      {
        GMRSS_MPEbash::cGMRSS_MPEarg arg;
          if(arg.Parse(argc, argv))
          {
            GMRSS_MPEbash::bash_run(argc, argv);
          }else
          {
            // needed to ensure appropriate OpenGL context is created for VTK rendering.
            QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());   // must be here
              // QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
            // QT Stuff
            QApplication app( argc, argv );

            // {
            //     // Install TouchBarProvider as application delegate
            //     TouchBarProvider *touchBarProvider = [[TouchBarProvider alloc] init];
            //     [touchBarProvider installAsDelegateForApplication:[NSApplication sharedApplication]];
            // }

            QApplication::setStyle("fusion");

            qInitResources_icons();

            MainWindow myMainWindow;
            //  myMainWindow.showMaximized();
            myMainWindow.show();


            // {
            //     // Install TouchBarProvider as window delegate
            //     NSView *view = reinterpret_cast<NSView *>(textEdit.winId());
            //     TouchBarProvider *touchBarProvider = [[TouchBarProvider alloc] init];
            //     [touchBarProvider installAsDelegateForWindow:view.window];
            // }
            return app.exec();
          }
      }else
      {
        GMRSS_MPEbash::bash_run(argc, argv);
      }
  #endif
 
 return 0;
}
