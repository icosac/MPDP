#ifndef MPMD_RS_HH
#define MPMD_RS_HH

// Library includes
#include <curve.hh>
#include <utils.hh>

#define DEFAULT_KMAX 1

// System includes
#include <exception>

#ifndef SAVE_DEBUG
#define SAVE_DEBUG(debug, var) if (debug != nullptr) { debug->push_back(var); }
#endif

struct RSSegment {
  double x, y;
  double thi, thf;
  double l, k;
  int dir;

  RSSegment() {}

  RSSegment(double _x, double _y, double _thi, double _thf, double _l, double _k, int _dir) :
      x(_x), y(_y), l(_l), thi(_thi), thf(_thf), k(_k), dir(_dir) {}

  std::string toString() const {
    std::stringstream ss;
    ss << "x: " << x << ", y: " << y << ", thi: " << thi << ", thf: " << thf << ", l: " << l << ", k: " << k << ", dir: " << dir;
    return ss.str();
  }

  friend std::ostream & operator << (std::ostream & os, const RSSegment & seg) {
    os << seg.toString();
    return os;
  }
};

class RS : public Curve {
private:
  LEN_T _L1 = 0.0, _L2 = 0.0, _L3 = 0.0, _L = 1e100;
  int _Nman = 0;
  int _Nseg = 0;

  K_T _kmax = 0.0; ///<The maximum curvature and the curvature for each part of the Reeds-Sheep path.
  std::vector<double> K = { 0,0,0,0,0 };
  std::vector<double> L = { 0,0,0,0,0 };
  std::vector<double> D = { 0,0,0,0,0 };
  std::vector<double> X = { 0,0,0,0,0,0 };
  std::vector<double> Y = { 0,0,0,0,0,0 };
  std::vector<double> TH = { 0,0,0,0,0,0 };
  std::vector<std::string> ManType = { "","","","","" };

  Configuration2 circleLine(double s, double dir, double kur, double kmax, Configuration2 c, int seg);

public:
  RS() :
      Curve(CURVE_TYPE::RS),
      _kmax(0) {}

  /*!
   * @brief Constructor or the Reed-Shepp curve.
   * @detail The constructor depends on the number of `params`:
   *  - 0: then `Kmax` is set to `DEFAULT_KMAX` and the number of the manoeuvre to -1 so that all the possibilities are
   *  considered and computed
   *  - 1: then `Kmax` is set to `params[0]` and the number of the manoeuvre to -1 so that all the possibilities are
   *  considered and computed
   *  - 2: then `Kmax` is set to `params[0]` and the number of the manoeuvre to `params[1]` so that only that particular
   *  manoeuvre is considered and computed
   *  - 5: then `Kmax` is set to `params[0]` and the number of the manoeuvre to `params[1]` so that only that particular
   *  manoeuvre is considered, plus `_L1`, `_L2` and `_L3` are set so that only the `buildRS` function needs to be
   *  called.
   * @param ci Initial configuration
   * @param cf Final configuration
   * @param params A vector of possible parameters for the curve:
   *  - [0] is the curvature `Kmax`
   *  - [1] is the number of the manoeuvre, default is -1 to try all the manoeuvres
   *  - [2] is `_L1`
   *  - [3] is `_L2`
   *  - [4] is `_L3`
   */
  RS(Configuration2 ci, Configuration2 cf, std::vector<real_type> params = {}) :
      Curve(ci, cf, CURVE_TYPE::RS, params)
  {
    if (params.empty()) {
      this->_kmax = DEFAULT_KMAX;
      this->_Nman = -1;
    }
    else if (params.size()==1) {
      this->_kmax = params[0];
      this->_Nman = -1;
    }
    else if (params.size()==2) {
      this->_kmax = params[0];
      this->_Nman = (int)params[1];
    }
    else {
      this->_kmax = params[0];
      this->_Nman = (int)params[1];
      this->_L1 = params[2];
      this->_L2 = params[3];
      this->_L3 = params[4];
    }
  }

  /**
   * @brief Function that computes the best manoeuvre for the given initial and final configuration.
   */
  void solve () {
    this->reeds_shepp(this->_Nman);
    this->buildRS(this->_Nman);
  }

  void buildRS(int man = -1);
  double reeds_sheppa(std::vector<double>* debug = nullptr);
  double reeds_shepp(int Nman = -1, std::vector<double>* debug = nullptr);

  [[nodiscard]] int getNseg()                         const { return this->_Nseg; }
  [[nodiscard]] int getNman()                         const { return this->_Nman; }
  [[nodiscard]] K_T getKmax()                         const { return this->_kmax; }
  [[nodiscard]] LEN_T getL1()                         const { return this->_L1; }
  [[nodiscard]] LEN_T getL2()                         const { return this->_L2; }
  [[nodiscard]] LEN_T getL3()                         const { return this->_L3; }
  [[nodiscard]] std::vector<double> getK()            const { return this->K; }
  [[nodiscard]] std::vector<double> getX()            const { return this->X; }
  [[nodiscard]] std::vector<double> getY()            const { return this->Y; }
  [[nodiscard]] std::vector<double> getL()            const { return this->L; }
  [[nodiscard]] std::vector<double> getTH()           const { return this->TH; }
  [[nodiscard]] std::vector<double> getD()            const { return this->D; }
  [[nodiscard]] std::vector<std::string> getManType() const { return this->ManType; }
  [[nodiscard]] LEN_T l() const override { return this->_L; }
  [[nodiscard]] LEN_T sumL() const
  {
    double l = 0.0;
    for (auto tmp : this->L){
      l+= std::abs(tmp);
    }
    return l;
  }

  [[nodiscard]] std::vector<RSSegment> getSegmentsData();

  [[nodiscard]] inline std::string getManTypeS() const {
    std::string ret = this->ManType[0];
    for(int i=1; i<this->getNseg(); i++){
      ret+=this->ManType[i];
    }
    return ret;
  }

  void setKmax (K_T k) { this->_kmax = k; }

  std::string to_string(int prec = 16) const
  {
    std::stringstream out;
    out << std::endl;
    out << std::setprecision(prec);
    out << "Man      = " << this->getNman() << std::endl;
    out << "Kmax     = " << this->getKmax() << std::endl;
    out << "Segments = " << this->getNseg() << std::endl;
    out << "tr       = " << this->getL1() << std::endl;
    out << "ur       = " << this->getL2() << std::endl;
    out << "vr       = " << this->getL3() << std::endl;
    out << "-------------------------------" << std::endl;
    out << "Sum len  = " << this->sumL() << std::endl;
    out << "Length   = " << this->l() << std::endl;
    out << "-------------------------------" << std::endl;
    for (int i = 0; i < this->getNseg(); ++i) {
      out << i+1 <<
          " (x,y,th) => " << this->getX()[i] << "  " << this->getY()[i] << "  " << this->getTH()[i] <<
          " to (x,y,th) => " << this->getX()[i+1] << "  " << this->getY()[i+1] << "  " << this->getTH()[i+1] <<
          " of type:  " << this->getManType()[i] << std::endl;
    }
    return out.str();
  }

  friend std::ostream& operator<< (std::ostream& out, const RS & curve)
  {
    out << curve.to_string();
    return out;
  }

#ifdef MPDP_DRAW
  void draw();
#endif
};

#endif //MPMD_RS_HH
