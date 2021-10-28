#include "asyplot.hh"



  AsyPlot::AsyPlot( string _filename )
  : filename(_filename)
  {
	if (!openFile()) throw std::runtime_error("Failed to open file " + filename);
    initFile();
  }

  AsyPlot::~AsyPlot() {
    if ( showAxes ) displayAxes();
    if ( closeFile() ) compileFile();
  }

  void
  AsyPlot::compileFile() {
    string cmdComp = "asy -f pdf "+ filename;
    system(cmdComp.c_str());
    string pdfFile = filename.substr(0,filename.find(".asy")) + ".pdf";
    std::cout << pdfFile << std::endl;
    string cmdOpen = "(okular " + pdfFile + " &> /dev/null )&";
    //system(cmdOpen.c_str());
  }

  void
  AsyPlot::initFile() {
    file
      << "// File generated automatically from C++ \n\n\n"
      << "import graph; import palette;\n"
      << "include \"clothoidLib.asylib\";\n"
      << "size(10cm,10cm);\n"
      << "\n\n\n" ;
  }

  
  void AsyPlot::writeLine(string const & line) const {
	  file << line << "\n";
  }

 


  void
  AsyPlot::dot( double x, double y, string const & penna ) const {
    file << "dot((" << x << "," << y << ")," << penna << ");\n\n";
  }

  
  void
  AsyPlot::drawRect( double x0, double y0,
                     double x1, double y1,
                     double x2, double y2,
                     double x3, double y3,
                     string const & penna ) const {
	file
    << "draw((" << x0 << "," << y0 << ") -- "
    << '(' << x1 << ',' << y1 << ") -- "
    << '(' << x2 << ',' << y2 << ") -- "
    << '(' << x3 << ',' << y3 << ") -- cycle, "
    << penna << ");\n\n";
  }

  void
  AsyPlot::fillRect(double x0, double y0,
	           	    double x1, double y1,
		            double x2, double y2,
					double x3, double y3,
					string const & pennaBordo,
					string const & pennaRiemp) const {
	  file
		  << "filldraw((" << x0 << "," << y0 << ") -- "
		  << '(' << x1 << ',' << y1 << ") -- "
		  << '(' << x2 << ',' << y2 << ") -- "
		  << '(' << x3 << ',' << y3 << ") -- cycle, "
		  << pennaRiemp << ", " << pennaBordo << "); \n\n";
  }

  void
  AsyPlot::displayAxes() const {
    file
      << "xaxis(\"$x$\", black+fontsize(7pt),xmin=0,xmax=5,Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n"
      << "yaxis(\"$y$\", black+fontsize(7pt),ymin=0,ymax=5,Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n";
  }

  void
  AsyPlot::displayAxes( string const & labX,
                        string const & labY,
                        double xmin,
                        double xmax,
                        double ymin,
                        double ymax ) const {
// 	  file
//       << "xaxis(\"" << labX << "\", black+fontsize(7pt),xmin="
//       << xmin << ",xmax=" << xmax
//       << ",Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n"
//       << "yaxis(\"" << labY << "\", black+fontsize(7pt),ymin="
//       << ymin << ",ymax=" << ymax
//       << ",Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n";

      file
      << "xaxis(\"" << labX << "\", BottomTop(), Ticks(Label(\"$%.2f$\"), Step=0.4,step=0.1, pTick=.8red, ptick=lightgrey, extend=true));\n"
      << "yaxis(\"" << labY << "\", LeftRight(), RightTicks(Label(\"$%.2f$\"), Step=0.4,step=0.1, pTick=.8red, ptick=lightgrey, extend=true));\n";
//       xaxis( "x [m]", BottomTop(), Ticks(Label("$%.2f$"), Step=0.8, step=0.1, pTick=.8red, ptick=lightgrey, extend=true));
// yaxis( "y [m]" , LeftRight(), RightTicks(Label("$%.2f$"), Step=0.8, step=0.1, pTick=.8red, ptick=lightgrey, extend=true));
  }

  void
  AsyPlot::xyLimits(double      xmin,
		            double      xmax,
	                double      ymin,
					double      ymax,
					bool           crop) const {
	  if (crop)
		  file << "xlimits(" << xmin << ", " << xmax << ", Crop);\n"
		       << "ylimits(" << ymin << ", " << ymax << ", Crop);\n";
	  else
		  file << "xlimits(" << xmin << ", " << xmax << ");\n"
		       << "ylimits(" << ymin << ", " << ymax << ");\n";
  }
  void
  AsyPlot::drawLine( double x0, double y0,
                     double x1, double y1,
                     string const & penna ) const {
    file
      << "draw((" << x0 << "," << y0 << ") -- "
      << "(" << x1 << "," << y1 << "), " << penna << ");\n\n";
  }

  void
  AsyPlot::drawAngle(double s0, double th0,
		             double s1, double th1, double dk,
		             string const & penna, int npts) const {
	  file << "path parabola = ";
	  double ds = (s1 - s0) / npts;
	  double s;
	  for (int i = 0; i < npts-1; ++i) {
		  s = s0 + i * ds;
		  file << "(" << s << ", " <<  ((s0 - s1)*(-s1 + s)*(s - s0)*dk + 2 * s0*th1 - 2 * s1*th0 + 2 * s*(th0 - th1)) / (2 * s0 - 2 * s1) << ")--";

	  }
	  s = s1;
	  file << "(" << s << ", " << ((s0 - s1)*(-s1 + s)*(s - s0)*dk + 2 * s0*th1 - 2 * s1*th0 + 2 * s*(th0 - th1)) / (2 * s0 - 2 * s1) << ");\n\n";
	  file << "draw(parabola, " << penna << ");\n\n";
  }

  void
  AsyPlot::label( string const & text,
                  double      x,
                  double      y,
                  string const & placement,
                  string const & penna ) const {
  	file
      << "label(\"" << text << "\", (" << x << ", " << y << "), "
  	  << (!placement.empty() ? placement + ", " : "")
	    << penna << ");\n\n";
  }

  void
  AsyPlot::draw(string const & text,
		  string const & penna) const {
	  file
		  << "draw(" << text << ", "
		  << penna << ");\n\n";
  }

  bool
  AsyPlot::openFile() {
    file.open(filename);
    return file.is_open();
  }

  bool
  AsyPlot::closeFile() {
    file.close();
    return !file.is_open();
  }


