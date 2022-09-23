#ifdef CUDA_ON
#include <dp.cuh>

namespace DP {
  namespace {
    class Cell {
    private:
      Angle _th;   ///<Angle of the final point of the point to point curve.
      LEN_T _l;    ///<Length of the point to point curve.
      int _nextID; ///<Id to the next cell for dynamic programming.

    public:
      /*!
       * Default void constructor which returns a cell initialized with ANGLE::FREE, max length and -1 as next cell.
       */
      Cell() : _th(ANGLE::FREE), _l(MAX_LEN_T), _nextID(0) {}

      /*!
       * Constructor that takes in input an angle, a length and the next id and returns a DP::Cell.
       * @param th The initial angle of the point to point curve.
       * @param l The length of the point to point curve.
       * @param next The next id of the cell.
       */
      BOTH Cell(Angle th, LEN_T l, int nextID) :
          _th(th), _l(l), _nextID(nextID) {}

      BOTH Angle th() const { return this->_th; }      ///<Returns the angle.
      BOTH LEN_T l()  const { return this->_l; }       ///<Returns the length.
      BOTH int next() const { return this->_nextID; }  ///<Returns the next id.

      /*!
       * Sets the new angle.
       * @param th The new angle to be set.
       * @return the new set angle.
       */
      BOTH Angle th(Angle th) {
        this->_th = th;
        return this->th();
      }
      /*!
       * Sets the new length.
       * @param th The new length to be set.
       * @return the new set length.
       */
      BOTH LEN_T l(LEN_T l) {
        this->_l = l;
        return this->l();
      }
      /*!
       * Sets the new next id.
       * @param th The new next id to be set.
       * @return the new set next id.
       */
      BOTH int next(int nextID){
        //printf("nextID in class %d %u\n", nextID, nextID);
        this->_nextID = nextID;
        return this->next();
      }

      /*!
       * Creates a deep copy of a cell to `this`.
       * @param d The cell to copy from.
       * @return `*this`.
       */
      BOTH Cell& copy(const Cell &d) {
        this->th(d.th());
        this->l(d.l());
        this->next(d.next());

        return *this;
      }
      /*!
       * Overrides the assign operator (=) to make a deep copy of a cell to `this`.
       * @param d The cell to copy from.
       * @return `*this`.
       */
      BOTH Cell& operator=(const Cell &d) {
        this->copy(d);
        return *this;
      }

      /*!
       * Function to print the most essential info about `DP::Cell`.
       * @param pretty An additional truth value to print a prettier version. Default is `false`.
       * @return A `std::stringstream` object containing the data of `DP::Cell`.
       */
      std::stringstream to_string(bool pretty = false) const {
        std::stringstream out;
        out << std::setw(20) << std::setprecision(17);
        if (pretty) {
          out << "th: " << this->th() << " l: " << this->l();
        } else {
          out << "<" << (Angle)(this->th()*1.0) << ", " << (LEN_T)(this->l()) << ">";
        }
        return out;
      }
      /*! This function overrides the << operator so to print with `std::cout` the most essential info about the `DP::Cell`.
          \param[in] out The out stream.
          \param[in] data The `DP::Cell` to print.
          \returns An output stream to be printed.
      */
      friend std::ostream &operator<<(std::ostream &out, const Cell &data) {
        out << data.to_string().str();
        return out;
      }

    };
  } //Anonymous namespace to hide information
} //Namespace DP

__global__ void dubinsWrapper(Configuration2 c0, Configuration2 c1, double Kmax, double* L){
  CURVE c(c0, c1, Kmax);
  //printf("%.17f\n", c.l());
  L[0]+=c.l();
}

__global__ void printResults(real_type* results, uint discr, uint size){
  for (int i=0; i<size; i++){
    for(int j=0; j<discr; j++){
      for(int h=0; h<discr; h++){
        printf("(%2.0f,%.2f)", (float)((i*discr+j)*discr+h), results[(i*discr+j)*discr+h]);
      }
      printf("\t");
    }
    printf("\n");
  }
}

__global__ void printMatrix(DP::Cell* matrix, uint discr, uint size){
  for (int i=0; i<size; i++){
    for(int j=0; j<discr; j++){
      printf("(%d,%d)", (i*discr+j), matrix[i*discr+j].next());
    }
    printf("\n");
  }
}


// returns (up to) two circles through two points, given the radius
// Marco Frego and Paolo Bevilacqua in "An Iterative Dynamic Programming Approach to the Multipoint Markov-Dubins Problem" 2020.
static inline
void circles(real_type x1, real_type y1, real_type x2, real_type y2, real_type r, std::vector<real_type> & XC, std::vector<real_type> & YC) 
{
  real_type TOL = 1e-8;
  
  real_type q = std::hypot(x2-x1, y2-y1);
  real_type x3 = 0.5*(x1+x2);
  real_type y3 = 0.5*(y1+y2);

  real_type delta = r*r-q*q/4.;
    
  XC.clear();
  YC.clear();

  if (delta < -TOL) {
    return;
  }
  
  if (delta < TOL) 
  {
    XC.push_back(x3);
    YC.push_back(y3);
  }
  else
  {
    real_type deltaS = std::sqrt(delta);
    XC.push_back(x3 + deltaS*(y1-y2)/q);
    YC.push_back(y3 + deltaS*(x2-x1)/q);
    XC.push_back(x3 - deltaS*(y1-y2)/q);
    YC.push_back(y3 - deltaS*(x2-x1)/q);
  }
}

// Marco Frego and Paolo Bevilacqua in "An Iterative Dynamic Programming Approach to the Multipoint Markov-Dubins Problem" 2020.
// The function name is pretty self-explainatory 
uint guessInitialAngles(std::vector<std::set<Angle> >& moreAngles, const std::vector<Configuration2>& points, const std::vector<bool> fixedAngles, const real_type K){
  uint max=0;
  for (uint i=1; i<points.size(); i++){
    moreAngles.push_back(std::set<Angle>());
    if (i==1) { moreAngles.push_back(std::set<Angle>()); }
    //First add the lines connecting two points:
    Angle th = std::atan2((points[i].y()-points[i-1].y()), (points[i].x()-points[i-1].x()));
    if (!fixedAngles[i-1]){ moreAngles[i-1].insert(th); }
    if (!fixedAngles[i])  { moreAngles[i].insert(th); }
    
    //Then add the possible angles of the tangents to two possible circles:
    std::vector<real_type> XC, YC;
    circles(points[i-1].x(), points[i-1].y(), points[i].x(), points[i].y(), 1./K, XC, YC);
    
    for (uint j=0; j<XC.size(); j++){
      if (!fixedAngles[i-1]){
        th = std::atan2(points[i-1].y()-YC[j], points[i-1].x()-XC[j]);
        moreAngles[i-1].insert(th+M_PI/2.);
        moreAngles[i-1].insert(th-M_PI/2.);
      }
      if (!fixedAngles[i]){
        th = std::atan2(points[i].y()-YC[j], points[i].x()-XC[j]);
        moreAngles[i].insert(th+M_PI/2.);
        moreAngles[i].insert(th-M_PI/2.);
      }
    }
    if (moreAngles[i-1].size()>max){
      max=moreAngles[i-1].size();
    }
    if (i==points.size()-1 && moreAngles[i].size()>max){
      max=moreAngles[i].size();
    }
  }  
  return max;
}

std::vector<Angle> bestAngles(DP::Cell* matrix, int discr, int size){
  DP::Cell* best=&matrix[0];
  //Find best path
  for (int i=size; i<discr*size; i+=size){
    if (best->l()>matrix[i].l()  && matrix[i].l()!=0){ //TODO The second check is actually a bug in solveCell, but I'm not in the right mind to find this bug, please fix later
      best=&matrix[i];
    }
  }
  //Retrieve best angles
  std::vector<Angle> ret(1, best->th());
  uint nextID=best->next();
  while (nextID!=0){
    ret.push_back(matrix[nextID].th());
    nextID=matrix[nextID].next();
  }
  return ret;
}

std::pair<LEN_T, std::vector<Angle> >
bestAnglesMatrix(DP::Cell* matrix, int discr, int size, const std::vector<bool>& fixedAngles){
  DP::Cell* best=&matrix[0];

  if (!fixedAngles[0]){
    for(int i=1; i<discr; i++){
      if (matrix[i].l()<best->l())
        best=&matrix[i];
    }
  }

  //std::cout << "In function Length: " << std::setw(20) << std::setprecision(17) << best->l() << std::endl;

  std::vector<real_type> ret={best->th()};
//  std::vector<real_type> ret={best->l(), best->th()};
  int nextID=best->next()+discr;
  for (int i=1; i<size; i++){
    ret.push_back(matrix[nextID].th());
    nextID=matrix[nextID].next()+(i+1)*discr;
  }
  //ret.insert(ret.begin(), best->l());
  return std::pair<LEN_T, std::vector<Angle> >(best->l(), ret);
}

__global__ void solveCol( DP::Cell* matrix, uint discr, uint size, const bool* fixedAngles, 
                          Configuration2 c0, Configuration2 c1, 
                          Angle a00, Angle a01, real_type* params, int i, Angle fullAngle, bool halveDiscr
                        ){
  int tidx=threadIdx.x+blockDim.x*blockIdx.x;
  int stride=blockDim.x*gridDim.x;
  int halfDiscr=(discr-1)/2;
  int j=tidx;

  // if (j<discr){
  for (; j<discr; j+=stride){
    Angle bestA=0.0;
    LEN_T bestL=MAX_LEN_T; 
    int bestK=0;
    if (!fixedAngles[i-1]){ //If angle is fixed I don't have to change it
      double hj=fullAngle*((j-halfDiscr)*1.0)/(((halveDiscr ? halfDiscr : discr)*1.0));
      c0.th(a00+hj); 
    } 
    
    for (int k=0; k<discr; k++){ //SolveCell
      LEN_T currL=MAX_LEN_T;
      if (!fixedAngles[i]){ //If angle is fixed I don't have to change its
        double hk=fullAngle*((k-halfDiscr)*1.0)/(((halveDiscr ? halfDiscr : discr)*1.0));
        c1.th(a01+hk); 
      } 
      CURVE c=CURVE(c0, c1, params); 
      DP::Cell* next=(i==size-1 ? NULL : &matrix[k*size+(i+1)]);
      if (c.l()>0){
        currL=c.l();
        if (next!=NULL){
          currL+=next->l();
        }  
        if (currL<bestL || bestL==MAX_LEN_T){
          bestL=currL;
          bestA=c1.th();
          bestK=k;
        }
      }
      if (fixedAngles[i]){ k=discr; } //If the angle is fixed I don't have to change it
    }
    
    if (bestL!=MAX_LEN_T){
      uint nextID=(i==size-1 ? 0 : bestK*size+(i+1));
      matrix[j*size+i]=DP::Cell(bestA, bestL, nextID);
    }
    if (i==1){
      matrix[size*j]=DP::Cell(c0.th(), bestL, (size*j+i));
    }
    if(fixedAngles[i-1]) j=discr;
  }
}


__global__ void solveMatrixCol (DP::Cell* matrix, uint discr, uint size, const bool* fixedAngles, 
                                Configuration2 c0, Configuration2 c1, 
                                real_type* params, int i, uint ref=0){
  uint tidx=threadIdx.x+blockDim.x*blockIdx.x;
  uint stride=blockDim.x*gridDim.x;

  uint j=tidx;
  // if (j<discr){
  for (; j<discr; j+=stride){
    c0.th(matrix[i*discr+j].th());
    for (int h=0; h<(int)(discr); h++){
      c1.th(matrix[(i+1)*discr+h].th());

      CURVE c=CURVE(c0, c1, params);
      LEN_T currL=c.l()+matrix[(i+1)*discr+h].l();
      //if (ref==3 && i==0 && j==0){
      //  printf("x0: %.2f y0: %.2f th0: %.16f x1: %.2f y1: %.2f th1: %.16f matrix[i*discr+j].l(): %.16f currL %.16f c.l(): %.16f matrix[(i+1)*discr+h].l(): %.16f\n", c0.x(), c0.y(), c0.th(), c1.x(), c1.y(), c1.th(), (matrix[i*discr+j].l()<10000.0 ? matrix[i*discr+j].l() : 10000.0), currL, c.l(), matrix[(i+1)*discr+h].l());
      //}
      if (currL<matrix[i*discr+j].l()) {
        matrix[i*discr+j].l(currL);
        //printf("nextID in func: %u %d\n", h, h);
        matrix[i*discr+j].next(h);
      }
      if (fixedAngles[i+1]) {h=discr;}
    }
    if (matrix[i*discr+j].next()==-1) printf("[%u] BIG NO\n", i*discr+j);
    if (fixedAngles[i]) {j=discr;}
  }
}

void solveDPMatrix (std::vector<Configuration2> points, DP::Cell* dev_matrix, uint discr, std::vector<bool> fixedAngles, 
                    bool* dev_fixedAngles, real_type* dev_params, uint nThreads=128, uint ref=0){

  //REMOVE
  size_t size=points.size();
  int numberOfSMs; cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, cudaGetdeviceID());


  for (int i=size-2; i>=0; i--){
    Configuration2 c0=points[i];
    Configuration2 c1=points[i+1];

    size_t threads=discr>nThreads ? nThreads : discr;
    size_t blocks=((int)(discr/threads)+1)*numberOfSMs; 
    
    if(fixedAngles[i]){
      threads=1;
      blocks=1;
    }
    solveMatrixCol<<<blocks, threads>>>(dev_matrix, discr, size, dev_fixedAngles, c0, c1, dev_params, i, ref);
    cudaDeviceSynchronize();
    cudaCheckError(cudaGetLastError()); 
  }
  if (ref==30){
    printMatrix<<<1,1>>>(dev_matrix, discr, size);
    cudaDeviceSynchronize();
  }
}

std::pair<LEN_T, std::vector<Angle> >
solveDPMatrixAllocator (std::vector<Configuration2> points, uint discr,
                        const std::vector<bool> fixedAngles, std::vector<real_type> params,
                        Angle fullAngle=2*M_PI, uint nThreads=0, uint ref=0){
  size_t size=points.size();
  DP::Cell* matrix;
  bool* dev_fixedAngles=cudaSTDVectorToArray<bool>(fixedAngles);
  real_type* dev_params=cudaSTDVectorToArray<real_type>(params);  
  DP::Cell* dev_matrix;

  int numberOfSMs; cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, cudaGetdeviceID());
  
  std::vector<std::set<Angle> > moreAngles;
  uint addedAngles=0;

  addedAngles=guessInitialAngles(moreAngles, points, fixedAngles, params[0]);

  uint halfDiscr=(uint)((discr-(discr%2==0 ? 0 : 1))/2);
  real_type dtheta=fullAngle/(((int)(discr/2))*1.0);
  if(ref==0){
    dtheta=fullAngle/(discr*1.0);
  }
  discr+=addedAngles;
  cudaMallocHost(&matrix, (sizeof(DP::Cell)*size*discr));
  cudaMalloc(&dev_matrix, sizeof(DP::Cell)*size*discr);

  for (uint i=0; i<size; i++){ //TODO change back, remove l=-1 if fixedAngles
    LEN_T l = (i==size-1 ? 0 : std::numeric_limits<LEN_T>::max());
    for (uint j=0; j<=halfDiscr; j++){
      if (fixedAngles[i]){
        matrix[i*discr+j]=DP::Cell(points[i].th(), l, -1);
        break;
      }
      else {
        if(j==0) { 
          matrix[i*discr+j]=DP::Cell(points[i].th(), l, -1); 
        }
        else{
          matrix[i*discr+j]          =DP::Cell(mod2pi(points[i].th()-(j*1.0)*dtheta), l, -1);
          matrix[i*discr+j+halfDiscr]=DP::Cell(mod2pi(points[i].th()+(j*1.0)*dtheta), l, -1); 
        }
      }
    }
    if (true){
      uint j=discr-addedAngles;
      if (!fixedAngles[i]){
        for (std::set<Angle>::iterator it=moreAngles[i].begin(); it!=moreAngles[i].end(); ++it){
          matrix[i*discr+j]=DP::Cell(*it, l, -1);
          j++;
        }
        for (; j<discr; j++){
          matrix[i*discr+j]=DP::Cell(0, l, -1);
        }
      }
    }
  }
    
  cudaMemcpy(dev_matrix, matrix, sizeof(DP::Cell)*size*discr, cudaMemcpyHostToDevice);
  cudaCheckError(cudaGetLastError());

  solveDPMatrix(points, dev_matrix, discr, fixedAngles, dev_fixedAngles, dev_params, nThreads, ref);

  cudaMemcpy(matrix, dev_matrix, sizeof(DP::Cell)*size*discr, cudaMemcpyDeviceToHost);
  cudaCheckError(cudaGetLastError());


#ifdef DEBUG
  //Retrieve angles
  cout << "Computing best angles" << endl;
#endif
  std::pair<LEN_T, std::vector<Angle> > ret=bestAnglesMatrix(matrix, discr, size, fixedAngles);
#ifdef DEBUG
  std::vector<Angle> bestA=ret.second;
  printV(bestA)
#endif
  
#ifdef DEBUG
  LEN_T Length=0.0;
  for (unsigned int i=bestA.size()-1; i>0; i--){
    points[i].th(bestA[i]);
    points[i-1].th(bestA[i-1]);
    CURVE c(points[i-1], points[i], params.data());
    Length+=c.l();
  }
  cout << "\tMatrix length: " << setprecision(20) << Length << " " << setprecision(12) << (Length-7.46756219733842652175326293218) << endl;
  cout << "Printing for Matlab" << endl;
  cout << "X=[";
  for (unsigned int i=0; i<points.size(); i++){ cout << points[i].x() << (i!=points.size()-1 ? ", " : "];\n"); }
  cout << "Y=[";
  for (unsigned int i=0; i<points.size(); i++){ cout << points[i].y() << (i!=points.size()-1 ? ", " : "];\n"); }
  cout << "th=[";
  for (unsigned int i=0; i<bestA.size(); i++){ cout << bestA[i] << (i!=bestA.size()-1 ? ", " : "];\n"); }
  cout << "KMAX: " << params[0] << endl;
#endif

  cudaFreeHost(matrix);

  cudaFree(dev_matrix);
  cudaFree(dev_params);
  cudaFree(dev_fixedAngles);

  return ret;
}


__global__ void computeMore(DP::Cell* matrix, real_type* results, const bool* fixedAngles,
                            real_type* params, const Configuration2* points, 
                            size_t jump, size_t discr, size_t size, size_t iter){ //TODO It may be possible to remove cmp and use iter and i to compute the position
  uint tidx=threadIdx.x+blockDim.x*blockIdx.x;

  uint j=tidx;
  if (j<discr*jump*discr){        //j must be less than the number of rows (jump) times the number of inner cells per cell, times the number of cells per row
    uint cell=(int)(tidx/discr);  //The big cell
    uint inCell=tidx%discr;       //The small cell inside the big cell
    uint cmpId=(int)(cell/discr); //The row in this call
    uint pos=iter*jump+cmpId;     //The row w.r.t. the whole matrix

    if (pos<size-1){ 
      Configuration2 c0=points[pos];
      Configuration2 c1=points[pos+1];

      if (!fixedAngles[pos])   {c0.th(matrix[cell+iter*jump*discr].th());}
      if (!fixedAngles[pos+1]) {c1.th(matrix[inCell+(pos+1)*discr].th());}

      CURVE c=CURVE(c0, c1, params);
      if (c.l()>0){
        results[cell*discr+inCell]=c.l();
      }
    }
  }
}


void bestAnglesPerCell( DP::Cell* matrix, real_type* results, const std::vector<bool> fixedAngles, 
                        size_t size, size_t discr, size_t iter, size_t jump, size_t _threads, size_t numberOfSMs){
  //MIND THAT WE START FROM THE LAST ROW OF THE MATRIX
  int startRowIDM=iter*jump;        //Given #jump rows, this is the id of the first row in the jump group
  uint lastRowIDM=iter*jump+jump-1;  //Given #jump rows, this is the id of the last row in the jump group
  lastRowIDM=(lastRowIDM<size ? lastRowIDM : size-1);
  for (int i=lastRowIDM; i>=startRowIDM; i--){ //Cycle through the rows
    if (i==(int)(size-1)){continue;}    //If it's the last row, then skip it.
    uint startCellM=i*discr;              //Given #jump rows, this is the id of the cell in the row in the jump group I'm considering
    //std::cout << "startRowIDM: " << startRowIDM << std::endl;
    //std::cout << "lastRowIDM: " << lastRowIDM << std::endl;
    //std::cout << "startCellM: " << startCellM << std::endl;
    //std::cout << "i: " << i << std::endl;
    for (uint cellIDM=startCellM; cellIDM<startCellM+discr; cellIDM++){ //Cycle through all the cells in the row  
      for (uint h=0; h<discr; h++){        //Each cell in matrix corresponds to discr cells in results. Each h is in results is the same as the next row
        uint cellIDR=cellIDM*discr+h-startRowIDM*discr*discr;
        double currL=results[cellIDR]+matrix[(i+1)*discr+h].l();
        // int a=-2;
        if(currL<matrix[cellIDM].l()){
          matrix[cellIDM].l(currL);
          matrix[cellIDM].next(h);
          // a=matrix[cellIDM].next(h);
          //std::cout << "a: " << a << std::endl;
        }
        //if (cellIDM>29 && cellIDM<45){
        //  std::cout << "i: " << i << std::endl;
        //  std::cout << "cellIDM: " << cellIDM << std::endl;
        //  std::cout << "cellIDR: " << cellIDR << std::endl;
        //  std::cout << "results[cellIDR]: " << results[cellIDR] << std::endl;
        //  std::cout << "currL: " << currL << std::endl;
        //  std::cout << "matrix[cellIDM].l(): " << matrix[cellIDM].l() << std::endl;
        //  std::cout << "a: " << a << std::endl;
        //  std::cout << "matrix[cellIDM].n(): " << matrix[cellIDM].next() << std::endl;
        //}
        //std::cout << "matrix[(i+1)*discr+h].l(): " << matrix[(i+1)*discr+h].l() << std::endl;
        if (fixedAngles[i+1]){ h=discr; }
      } 
      if (matrix[cellIDM].next()<0) {printf("[%u] BIG NO\n", cellIDM);}
    }
    if (fixedAngles[i]){ i=startRowIDM-1; }
  }
}

std::pair<LEN_T, std::vector<Angle> >
solveDPAllIn1 ( std::vector<Configuration2> points, uint discr, const std::vector<bool> fixedAngles,
                std::vector<real_type> params, Angle fullAngle, uint nThreads=0, uint ref=0){

  //Get the number of multiproccessors in the GPU so to best compute the number of blocks afterwards.
  int numberOfSMs; cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, cudaGetdeviceID());


  uint addedAngles=0;
  std::vector<std::set<Angle> > moreAngles;
  addedAngles=guessInitialAngles(moreAngles, points, fixedAngles, params[0]);

  size_t size=points.size();
  //discr=(discr%2==0 ? discr+1 : discr); //So.... since we add always the angle in position 0, we'll always have an odd number of discretizionations... I'm not so sure about this, but ok
  uint halfDiscr=(uint)(discr/2);
  real_type dtheta=fullAngle/(((int)(discr/2))*1.0);
  if(ref==0){
    dtheta=fullAngle/(discr*1.0);
  }
  discr+=addedAngles;

  DP::Cell* matrix;
  cudaMallocHost(&matrix, sizeof(DP::Cell)*size*(discr));
  DP::Cell* dev_matrix;
  cudaMalloc(&dev_matrix, sizeof(DP::Cell)*size*(discr));

  bool* dev_fixedAngles=cudaSTDVectorToArray<bool>(fixedAngles);
  real_type* dev_params=cudaSTDVectorToArray<real_type>(params);  
  Configuration2* dev_points=cudaSTDVectorToArray<Configuration2>(points);

  for (uint i=0; i<size; i++){
    LEN_T l = (i==size-1 ? 0 : std::numeric_limits<LEN_T>::max());
    if (fixedAngles[i]){
      for (uint j=0; j<discr; j++){
        matrix[i*discr+j]=DP::Cell(points[i].th(), l, -1);
        //In this case I need to have the row full of the same values. Otherwise I should change the kernel function and add particular cases for fixed angles
      }
    }
    else{
      for (uint j=0; j<=halfDiscr; j++){
        COUT(j)
        if(j==0) { 
          matrix[i*discr+j]=DP::Cell(points[i].th(), l, -1); 
        }
        else{
          matrix[i*discr+j]          =DP::Cell(mod2pi(points[i].th()-(j*1.0)*dtheta), l, -1);
          matrix[i*discr+j+halfDiscr]=DP::Cell(mod2pi(points[i].th()+(j*1.0)*dtheta), l, -1); 
        }
      }
      uint j=discr-addedAngles;
      for (std::set<Angle>::iterator it=moreAngles[i].begin(); it!=moreAngles[i].end(); ++it){
        matrix[i*discr+j]=DP::Cell(*it, l, -1);
        j++;
      }
      for (; j<discr; j++){
        matrix[i*discr+j]=DP::Cell(points[i].th(), l, -1);
      }
    }
  }

  cudaMemcpy(dev_matrix, matrix, sizeof(DP::Cell)*size*discr, cudaMemcpyHostToDevice);
  cudaCheckError(cudaGetLastError());

  size_t jump=(params.size()>1 ? params[1] : 3);
  size_t iter=0;
  if ((size-1)%jump==0) { iter=(size-1)/jump; }
  else                  { iter=(size_t)(((size-1)+jump)/jump); }

  size_t totThreads=jump*discr*discr;
  size_t threads=totThreads>nThreads ? nThreads : totThreads;
  size_t blocks=((int)(totThreads/threads)+1)*numberOfSMs; 
  
  real_type *results, *dev_results1, *dev_results2, *dev_resultsapp;
  cudaMallocHost(&results, sizeof(real_type)*jump*discr*discr);
  cudaMalloc(&dev_results1, sizeof(real_type)*jump*discr*discr);
  cudaMalloc(&dev_results2, sizeof(real_type)*jump*discr*discr);
    
  for (int i=iter-1; i>=0; i--){
    computeMore<<<blocks, threads>>>(dev_matrix, dev_results1, dev_fixedAngles, dev_params, dev_points, jump, discr, size, i);
    cudaDeviceSynchronize();
    cudaCheckError(cudaGetLastError());
    
    cudaMemcpy(results, dev_results1, sizeof(real_type)*jump*discr*discr, cudaMemcpyDeviceToHost);
    cudaCheckError(cudaGetLastError());

    dev_resultsapp=dev_results1;
    dev_results1=dev_results2;
    dev_results2=dev_resultsapp;

    bestAnglesPerCell(matrix, results, fixedAngles, size, discr, i, jump, nThreads, numberOfSMs);
    cudaMemcpy(dev_matrix, matrix, sizeof(real_type)*size*discr, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    cudaCheckError(cudaGetLastError());
    
    #ifdef DEBUG
    printf("\n");
    #endif
  }
#ifdef DEBUG
  //Retrieve angles
  cout << "Computing best angles" << endl;
#endif
  std::pair<LEN_T, std::vector<Angle> > ret=bestAnglesMatrix(matrix, discr, size, fixedAngles);
#ifdef DEBUG
  std::vector<Angle> bestA=ret.second;
  printV(bestA)
#endif
  
#ifdef DEBUG
  LEN_T Length=0.0;
  for (unsigned int i=bestA.size()-1; i>0; i--){
    points[i].th(bestA[i]);
    points[i-1].th(bestA[i-1]);
    CURVE c(points[i-1], points[i], params.data());
    Length+=c.l();
  }
  cout << "\tAllInOne length: " << setprecision(20) << Length << " " << setprecision(12) << (Length-7.467562181965) << endl;

  cout << "Printing for Matlab" << endl;
  cout << "X=[";
  for (unsigned int i=0; i<points.size(); i++){ cout << points[i].x() << (i!=points.size()-1 ? ", " : "];\n"); }
  cout << "Y=[";
  for (unsigned int i=0; i<points.size(); i++){ cout << points[i].y() << (i!=points.size()-1 ? ", " : "];\n"); }
  cout << "th=[";
  for (unsigned int i=0; i<bestA.size(); i++){ cout << bestA[i] << (i!=bestA.size()-1 ? ", " : "];\n"); }
  cout << "KMAX: " << params[0] << endl;
#endif
  cudaFreeHost(matrix);
  cudaFreeHost(results);

  cudaFree(dev_matrix);
  cudaFree(dev_params);
  cudaFree(dev_fixedAngles);
  cudaFree(dev_points);
  cudaFree(dev_results1);
  cudaFree(dev_results2);

  return ret;
}

std::pair<LEN_T, std::vector<Angle> >
DP::solveDP(std::vector<Configuration2>& points, const std::vector<bool> fixedAngles,
            std::vector<real_type> params, int discr, uint nRefs, bool saveAngles,
            short type, uint threads, Angle _fullAngle){
  ////std::cout << "ciao2.1\n";
  if (points.size()!=fixedAngles.size()){
    std::cerr << "Number of points and number of fixed angles are not the same: " << points.size() << "!=" << fixedAngles.size() << std::endl;
    return std::pair<LEN_T, std::vector<Angle> >(MAX_LEN_T, std::vector<Angle>());
  }

  Angle fullAngle=_fullAngle;
  std::pair<LEN_T, std::vector<Angle> > ret;

  for(uint i=0; i<nRefs+1; ++i){
    switch(type){
      case 1:{
        ret=solveDPMatrixAllocator (points, discr, fixedAngles, params, fullAngle, threads, i);
        break;
      }
      case 2: default:{
        ret=solveDPAllIn1          (points, discr, fixedAngles, params, fullAngle, threads, i);
      }
    }
    std::vector<Angle> angles=ret.second;

    if(saveAngles){
      for (uint j=0; j<angles.size(); j++){
        if (!fixedAngles[j]){
          points[j].th(angles[j]);
        }
      }
    }

    if (i==0){
      fullAngle=fullAngle/(discr)*1.5;
      discr++; //This is because, yes.
    }
    else{
      fullAngle=fullAngle/(discr-1)*1.5;
    }
    angles.clear();
  }
  return ret;
}

#endif //CUDA_ON


