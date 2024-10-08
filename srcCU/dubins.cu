#ifdef CUDA_ON 

#include<dubins.cuh>

/*!
 * Function to scale to a standard settings the values. Credit to Marco Frego & Paolo Bevilacqua.
 */
BOTH void Dubins::scaleToStandard(Angle& phi, real_type& lambda, Angle& sth0, Angle& sth1, K_T& sKmax){
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
BOTH void Dubins::computeBest(Angle th0, Angle th1, real_type lambda, K_T& sKmax){
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

#ifdef __CUDA_ARCH__
  real_type len = MAX_LEN_T;
#else
  real_type len = std::numeric_limits<LEN_T>::max();
#endif 
  real_type temp1, temp2, temp3, t1, t2, t3, lc;

  // LSL
  real_type C = cos_1 - cos_0;
  real_type S = 2.0*sKmax + dsin;
  temp1 = atan2(C, S);
  temp2 = 2 + 4*Ksq - 2*dcos + 4*sKmax*dsin;
  if (temp2 >= 0) {
    temp3 = invK * sqrt(temp2);
    t1    = invK * mod2pi(temp1-th0);
    t2    = temp3;
    t3    = invK * mod2pi(th1-temp1);
    lc    = t1+t2+t3;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = 1; sk2 = 0; sk3 = 1;
      this->type(D_TYPE::LSL);
    }
  }

  // RSR
  C = -C;
  S = 2*sKmax - dsin;
  temp1 = atan2(C, S);
  temp2 = 2 + 4*Ksq - 2*dcos - 4*sKmax*dsin;
  if (temp2 >= 0) {
    temp3 = sqrt(invK*invK*(2 + 4*Ksq - 2*dcos - 4*sKmax*dsin));
    t1    = invK * mod2pi(th0-temp1);
    t2    = temp3;
    t3    = invK * mod2pi(temp1-th1);
    lc    = t1+t2+t3;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = -1; sk2 = 0; sk3 = -1;
      this->type(D_TYPE::RSR);
    }
  }

  // LSR
  C = scos;
  S = 2*sKmax + ssin;
  temp1 = atan2(-C, S);
  temp2 = -2 + 4*Ksq + 2*dcos + 4*sKmax*ssin;
  if (temp2 >= 0) {
    t2    = invK * sqrt(temp2);
    temp3 = -atan2(-2.0, (t2*sKmax));
    t1    = invK * mod2pi(-th0 + temp1 + temp3);
    t3    = invK * mod2pi(-th1 + temp1 + temp3);
    lc    = t1+t2+t3;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = 1; sk2 = 0; sk3 = -1;
      this->type(D_TYPE::LSR);
    }
  }

  // RSL
  // C = C
  S = 2*sKmax - ssin;
  temp1 = atan2(C, S);
  temp2 = -2 + 4*Ksq + 2*dcos - 4*sKmax*ssin;
  if (temp2 >= 0) {
    t2    = invK * sqrt(temp2);
    temp3 = atan2(2.0, (t2*sKmax));
    t1    = invK * mod2pi(th0 - temp1 + temp3);
    t3    = invK * mod2pi(th1 - temp1 + temp3);
    lc    = t1+t2+t3;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = -1; sk2 = 0; sk3 = 1;
      this->type(D_TYPE::RSL);
    }
  }

  // RLR
  C = dcos2;
  S = 2*sKmax - dsin;
  temp1 = atan2(C, S);
  temp2 = 0.125 * (6 - 4*Ksq  + 2*dcos + 4*sKmax*dsin);
  if (ABS<real_type>(temp2, 0.0) <= 1) {
    t2 = invK * mod2pi(2.0*M_PI - acos(temp2));
    t1 = invK * mod2pi(th0 - temp1 + 0.5*t2*sKmax);
    t3 = invK * mod2pi(dth+(t2-t1)*sKmax);
    lc = t1+t2+t3;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = -1; sk2 = 1; sk3 = -1;
      this->type(D_TYPE::RLR);
    }
  }

  // LRL
  C = -C;
  S = 2*sKmax + dsin;
  temp1 = atan2(C, S);
  temp2 = 0.125*(6 - 4*Ksq + 2*dcos - 4*sKmax*dsin);
  if (ABS<real_type>(temp2, 0.0) <= 1) {
    t2 = invK * mod2pi(2.0*M_PI - acos(temp2));
    t1 = invK * mod2pi(-th0 + temp1 + 0.5*t2*sKmax);
    t3 = invK * mod2pi(-dth + (t2-t1)*sKmax);
    lc = t1+t2+t3;
    if (lc < len) {
      len = lc; ss1 = t1; ss2 = t2; ss3 = t3;
      sk1 = 1; sk2 = -1; sk3 = 1;
      this->type(D_TYPE::LRL);
    }
  }

  //ScaleFromStandard
  this->s1(ss1*lambda);
  this->k1(sk1*this->kmax());
  this->s2(ss2*lambda);
  this->k2(sk2*this->kmax());
  this->s3(ss3*lambda);
  this->k3(sk3*this->kmax());
}

#endif //CUDA_ON