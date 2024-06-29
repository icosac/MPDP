#ifndef CUDA_ON
#include<dubins.hh>


static Configuration2 circleLine(double s, double kur, Configuration2 c) {
  double xEnd, yEnd, thetaEnd;
  xEnd = c.x() - f(s, kur , mod2pi(c.th() + m_pi));
  yEnd = c.y() - g(s, kur , mod2pi(c.th() + m_pi));
  thetaEnd = mod2pi(c.th() + kur  * s);

  Configuration2 res(xEnd, yEnd, thetaEnd, kur );

  return res;
}

/*!
 * Function to scale to a standard settings the values. Credit to Marco Frego & Paolo Bevilacqua.
 */
void Dubins::scaleToStandard(Angle& phi, real_type& lambda, Angle& sth0, Angle& sth1, K_T& sKmax){
  real_type dx = this->cf()->x() - this->ci()->x();
  real_type dy = this->cf()->y() - this->ci()->y();
  phi = atan2(dy, dx);
  lambda = hypot(dx, dy)*0.5;
  sKmax = this->kmax() * lambda;
  sth0 = mod2pi(this->ci()->th() - phi);
  sth1 = mod2pi(this->cf()->th() - phi);
}

/*!
 * Given the standardized version, compute the best word. Credit to Marco Frego & Paolo Bevilacqua.
 * @param th0 The initial standardized angle.
 * @param th1 The final standardized angle.
 * @param lambda A multiplier.
 * @param sKmax The standardized curvature.
 */
void Dubins::computeBest(Angle th0, Angle th1, real_type lambda, K_T& sKmax){
  K_T sk1=0.0, sk2=0.0, sk3=0.0;
  LEN_T ss1=0.0, ss2=0.0, ss3=0.0;

  real_type invK  = real_type(1)/sKmax;
  real_type sin_0 = sin(th0);
  real_type cos_0 = cos(th0);
  real_type sin_1 = sin(th1);
  real_type cos_1 = cos(th1);

  real_type Ksq   = sKmax*sKmax;
  real_type dcos  = cos(th0 - th1);
  real_type dcos2 = cos_0 - cos_1;
  real_type dsin  = sin_0 - sin_1;
  real_type scos  = cos_0 + cos_1;
  real_type ssin  = sin_0 + sin_1;

	real_type dth   = th0 - th1;

  real_type len = std::numeric_limits<real_type>::max();
  real_type temp1, temp2, temp3, t1, t2, t3, lc;

  // LSR
  real_type C = scos;
  real_type S = 2*sKmax + ssin;
  temp1 = std::atan2(-C, S);
  temp2 = -2 + 4*Ksq + 2*dcos + 4*sKmax*ssin;
  if (temp2 >= 0) {
    t2    = invK * std::sqrt(temp2);
    temp3 = -std::atan2(-2, (t2*sKmax));
    t1    = invK * mod2pi(-th0 + temp1 + temp3);
    t3    = invK * mod2pi(-th1 + temp1 + temp3);
    lc    = t1+t2+t3;
    std::cout << std::setprecision(12) << "LSR len: " << lc << std::endl;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = 1; sk2 = 0; sk3 = -1;
      this->dtype(D_TYPE::LSR);
    }
  }
  else {
    std::cout << "Cannot find valid LSR" << std::endl;
  }

  // LSL
  C = cos_1 - cos_0;
  S = 2*sKmax + dsin;
  temp1 = std::atan2(C, S);
  temp2 = 2 + 4*Ksq - 2*dcos + 4*sKmax*dsin;
  if (temp2 >= 0) {
    temp3 = invK * std::sqrt(temp2);
    t1    = invK * mod2pi(temp1-th0);
    t2    = temp3;
    t3    = invK * mod2pi(th1-temp1);
    lc    = t1+t2+t3;
    std::cout << std::setprecision(12) << "LSL len: " << lc << std::endl;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = 1; sk2 = 0; sk3 = 1;
      this->dtype(D_TYPE::LSL);
    }
  }
  else {
    std::cout << "Cannot find valid LSL" << std::endl;
  }

  // RSR
  C = -C;
  S = 2*sKmax - dsin;
  temp1 = std::atan2(C, S);
  temp2 = 2 + 4*Ksq - 2*dcos - 4*sKmax*dsin;
  if (temp2 >= 0) {
    temp3 = invK * std::sqrt(temp2);
    t1    = invK * mod2pi(th0-temp1);
    t2    = temp3;
    t3    = invK * mod2pi(temp1-th1);
    lc    = t1+t2+t3;
    std::cout << std::setprecision(12) << "RSR len: " << lc << std::endl;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = -1; sk2 = 0; sk3 = -1;
      this->dtype(D_TYPE::RSR);
    }
  }
  else {
    std::cout << "Cannot find valid RSR" << std::endl;
  }

  // RSL
  // C = C
  S = 2*sKmax - ssin;
  temp1 = std::atan2(C, S);
  temp2 = -2 + 4*Ksq + 2*dcos - 4*sKmax*ssin;
  if (temp2 >= 0) {
    t2    = invK * std::sqrt(temp2);
    temp3 = std::atan2(2, (t2*sKmax));
    t1    = invK * mod2pi(th0 - temp1 + temp3);
    t3    = invK * mod2pi(th1 - temp1 + temp3);
    lc    = t1+t2+t3;
    std::cout << std::setprecision(12) << "RSL len: " << lc << std::endl;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = -1; sk2 = 0; sk3 = 1;
      this->dtype(D_TYPE::RSL);
    }
  }
  else {
    std::cout << "Cannot find valid RSL" << std::endl;
  }

  // RLR
  C = dcos2;
  S = 2*sKmax - dsin;
  temp1 = std::atan2(C, S);
  temp2 = 0.125 * (6 - 4*Ksq  + 2*dcos + 4*sKmax*dsin);
  if (std::abs(temp2) <= 1) {
    t2 = invK * mod2pi(2*m_pi - std::acos(temp2));
    t1 = invK * mod2pi(th0 - temp1 + 0.5*t2*sKmax);
    t3 = invK * mod2pi(dth+(t2-t1)*sKmax);
    lc = t1+t2+t3;
    std::cout << std::setprecision(12) << "RLR len: " << lc << std::endl;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = -1; sk2 = 1; sk3 = -1;
      this->dtype(D_TYPE::RLR);
    }
  }
  else {
    std::cout << "Cannot find valid RLR" << std::endl;
  }

  // LRL
  C = -C;
  S = 2*sKmax + dsin;
  temp1 = std::atan2(C, S);
  temp2 = 0.125*(6 - 4*Ksq + 2*dcos - 4*sKmax*dsin);
  if (std::abs(temp2) <= 1) {
    t2 = invK * mod2pi(2*m_pi - std::acos(temp2));
    t1 = invK * mod2pi(-th0 + temp1 + 0.5*t2*sKmax);
    t3 = invK * mod2pi(-dth + (t2-t1)*sKmax);
    lc = t1+t2+t3;
    std::cout << std::setprecision(12) << "LRL len: " << lc << std::endl;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = 1; sk2 = -1; sk3 = 1;
      this->dtype(D_TYPE::LRL);
    }
  }
  else {
    std::cout << "Cannot find valid LRL" << std::endl;
  }

  //ScaleFromStandard
  this->s1(ss1*lambda);
  this->k1(sk1*this->kmax());
  this->s2(ss2*lambda);
  this->k2(sk2*this->kmax());
  this->s3(ss3*lambda);
  this->k3(sk3*this->kmax());
}

std::vector<std::vector<double>> Dubins::split_wise() {
  std::vector<std::vector<double>> res;
  uint n_split = 5;
  std::vector<double> tmp;

  Configuration2 curr = *this->ci();  
  for (uint seg=1; seg<4; seg++){
    // If the segment has length 0, skip it
    if (this->L(seg) == 0) continue;

    // If the segment is a straight line, add only the beginning and the end
    Configuration2 next_goal_c = circleLine(this->L(seg), this->k(seg), curr);
    if (this->k(seg) == 0 && next_goal_c.th() == curr.th()){
      tmp = {curr.x(), curr.y(), curr.th(), this->k(seg), this->L(seg)};
      res.push_back(tmp);
      curr = circleLine(this->L(seg), this->k(seg), curr);
    }
    // If the segment is a circle, add n_split points
    else {
      // Dynamically change the value of n_split so that the points are at least 0.01 apart
      n_split = std::min((uint) 3, (uint) std::ceil(this->L(seg)/0.01));

      for (uint i=0; i<n_split; i++){
        tmp = {curr.x(), curr.y(), curr.th(), this->k(seg), (this->L(seg)/n_split)};
        res.push_back(tmp);
        curr = circleLine(this->L(seg)/n_split, this->k(seg), curr);
      }
    }
  }
  
  return res;
}


#ifdef MPDP_DRAW
void Dubins::draw(std::ofstream& file, std::string label, size_t width, size_t height, bool solve, bool close, bool init) {
  if (solve) {
    this->solve();
  }

  if (init) {
    initAsyFile(file);
  }

  // Initial point
  Configuration2 c = this->ci()[0];
  file  << "p = clothoidPoints((" << c.x() << "," << c.y() << "), " << c.th()
        << "," << this->k1() << ", 0, " << this->s1() << ");" << std::endl;
  file << "draw(p,royalblue);" << std::endl;
  file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;
  if (label != ""){
    file << "label(\"$" << label << "$\", (" << c.x() << ", " << c.y() << "));" << std::endl;
  }

  // Intermediate point
  c = circleLine(this->s1(), this->k1(), c);
  file  << "p = clothoidPoints((" << c.x() << "," << c.y() << "), " << c.th()
        << "," << this->k2() << ", 0, " << this->s2() << ");" << std::endl;
  file << "draw(p,royalblue);" << std::endl;
  file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;

  // Final point
  c = circleLine(this->s2(), this->k2(), c);
  file  << "p = clothoidPoints((" << c.x() << "," << c.y() << "), " << c.th()
        << "," << this->k3() << ", 0, " << this->s3() << ");" << std::endl;
  file << "draw(p,royalblue);" << std::endl;
  file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;

  // Plot points
  file << "dot((" << this->ci()->x() << "," << this->ci()->y() << "), black);" << std::endl;
  file << "dot((" << this->cf()->x() << "," << this->cf()->y() << "), purple+3bp);" << std::endl;

  if (close){
    file.close();
  }
}
#endif //MPDP_DRAW


#endif //CUDA_ON