#ifndef Physics_Analysis_RscanDQ_H
#define Physics_Analysis_RscanDQ_H

class RscanDQ
{

public:

  RscanDQ(int runNo);
  ~RscanDQ();

public:

  int getStatus() {return m_status;}

  double getEbeam() {return m_Ebeam;}

  int getNumber() {return m_number;}

private:

  void Status(int runNo);

  void Ebeam(int runNo);

  void Number(int runNo);

private:

  int m_status;

  double m_Ebeam;

  int m_number;

};

#endif
