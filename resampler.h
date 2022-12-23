//
// Created by n.titov on 17.11.2022.
//

#ifndef LAGRANGE_CPP_RESAMPLER_H
#define LAGRANGE_CPP_RESAMPLER_H



class resampler {
public:
//  resampler();

  enum APPROXIMATION_TYPE{
    lagrange,
    sinq
  };
  virtual int do_approximation() = 0;
private:


};

//resampler* res = new resampler();

#endif // LAGRANGE_CPP_RESAMPLER_H
