#ifndef ParticleID_MRPCPID_H
#define ParticleID_MRPCPID_H
//
// MrpcPID package: particle Identification with Endcap MRPC detector
//


#include "ParticleID/ParticleIDBase.h"

class MrpcPID : public ParticleIDBase {

 public:
  static MrpcPID *instance();
  ~MrpcPID(){;} 

  void init();
  void calculate();
  bool IsPidInfoValid() const {return (m_ndof > 0); }
  double chi(int n) const {return m_chi[n];}
  double prob(int n) const {return m_prob[n];}
  double sigma(int n) const{return m_sigma[n];}
  double offset(int n) const{return m_offset[n];}
  int ndof() const {return m_ndof;}
  double mass2() const {return m_mass2;}
  int part() const {return m_part;}
  double rhit() const {return m_rhit;}
  int neuronPID() const {return -1;}
 protected:

  int neuronPIDCalculation() { return -1;}
  int particleIDCalculation();
  int LikelihoodCalculation() {return -1;}
  
 private:
  double m_chi[5];
  double m_prob[5];
  double m_sigma[5];
  double m_offset[5];
  double m_chimin;
  double m_pdfmin;
  int m_ndof;
  double m_mass2;
  int m_part; 
  double m_rhit; 
 private:
  MrpcPID();
  static MrpcPID *m_pointer;
};

#endif
