#ifndef MAGNETICFIELDSVC_H
#define MAGNETICFIELDSVC_H

// Include files
#include <vector>
#include <string>
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
//  backwards compatibility will be enabled ONLY in CLHEP 1.9
typedef HepGeom::Point3D<double> HepPoint3D;
typedef HepGeom::Vector3D<double> HepVector3D;
#endif
#include "MagneticField/ConnectionDB.h"

#ifndef BEAN
#include "GaudiKernel/Service.h"
//#include "GaudiKernel/IMagneticFieldSvc.h"
#include "MagneticField/IMagneticFieldSvc.h"
#include "GaudiKernel/IIncidentListener.h"
using namespace CLHEP;

// Forward declarations
template <class TYPE> class SvcFactory;

class IMagneticFieldSvc;
class MucMagneticField;
class IDataProviderSvc;
#else
class MucMagneticField;
typedef bool StatusCode;
#endif

/** @class MagneticFieldSvc MagneticFieldSvc.h
 *  A service for finding the magnetic field vector at a given
 *  point in space.
 *
 */

#ifndef BEAN
class MagneticFieldSvc : public extends<Service, IMagneticFieldSvc>,
                         virtual public IIncidentListener {
//protected:
public:  
  /// Standard Constructor.
  /// @param  name   String with service name
  /// @param  svc    Pointer to service locator interface
  MagneticFieldSvc( const std::string& name, ISvcLocator* svc );

  /// Virtual destructor.
  virtual ~MagneticFieldSvc();

public:

  /// Initialise the service (Inherited Service overrides)
  virtual StatusCode initialize();

  void init_params();
  void init_params(std::vector<double> current, std::vector<double> beamEnergy, int runNo);
  /// Finalise the service.
  virtual StatusCode finalize();

  /** Query the available interfaces.
   * @param riid Requested interface ID
   * @param ppvInterface Pointer to requested interface
   * @return StatusCode indicating SUCCESS or FAILURE.
   */
  //virtual StatusCode queryInterface( const InterfaceID& riid, 
  //                                   void** ppvInterface      );

  /// Service type.
  virtual const InterfaceID& type() const { return IMagneticFieldSvc::interfaceID(); };


  // Incident handler
  void handle(const Incident&);

#else
// -------------------------- BEAN ------------------------------------
class MagneticFieldSvc {

public:
  static MagneticFieldSvc* instance() {
    return (m_field) ? m_field : (m_field=new(std::nothrow)MagneticFieldSvc());
  }

private:
  MagneticFieldSvc();
  virtual ~MagneticFieldSvc();

public:
  void SetPath(std::string new_path)            {path = new_path;}
  std::string GetPath()                         {return path;}
  void SetPaths(std::string new_path, std::string new_DBFilePath)
      {this->SetPath(new_path); m_connect_run->SetDBFilePath(new_DBFilePath);}

  // function replacement for declareProperty()
  void TurnOffField (bool turnOffField = false) {m_turnOffField=turnOffField;}
  void GridDistance (int gridDistance = 5) {m_gridDistance = gridDistance;}
  void RunMode (int runmode = 2) {m_runmode = runmode;}
  void IfRealField (bool ifRealField = true) {m_ifRealField = ifRealField;}
  void OutLevel (int outlevel = 1) {m_outlevel = outlevel;}
  void SetScale (double scale = 1.0) {m_scale = scale;}
  void UniField (bool uniField = false) {m_uniField = uniField;}

  void Cur_SCQ1_55 (double Cur_SCQ1_55 = 349.4) {m_Cur_SCQ1_55 = Cur_SCQ1_55;}
  void Cur_SCQ1_89 (double Cur_SCQ1_89 = 426.2) {m_Cur_SCQ1_89 = Cur_SCQ1_89;}
  void Cur_SCQ2_10 (double Cur_SCQ2_10 = 474.2) {m_Cur_SCQ2_10 = Cur_SCQ2_10;}

  void UseDBFlag (bool useDB = true) {m_useDB = useDB;}
//   void ReadOneTime (bool readOneTime = false) {m_readOneTime = readOneTime;}
  void RunFrom (int runFrom = 8093) {m_runFrom = runFrom;}
  void RunTo (int runTo = 9025) {m_runTo = runTo;}

  /// Initialize
  bool init_mucMagneticField();
  virtual StatusCode initialize();
  void init_params(int runNo);

  // handler for new run
  void handle(int new_run);

  // get beam-energy vector
  const std::vector<double>& GetBeamEnergy() const {return beamEnergy;}
#endif

  /** IMagneticFieldSvc interface.
   * @param[in]  xyz Point at which magnetic field vector will be given
   * @param[out] fvec Magnectic field vector.
   * @return StatusCode SUCCESS if calculation was performed.
   */
  virtual StatusCode fieldVector( const HepPoint3D& xyz,
                                  HepVector3D& fvec ) const;

  virtual StatusCode uniFieldVector( const HepPoint3D& xyz,
                                  HepVector3D& fvec ) const;

  virtual double getReferField();

  virtual bool ifRealField() const;


private:
#ifndef BEAN
  /// Allow SvcFactory to instantiate the service.
  //friend class SvcFactory<MagneticFieldSvc>;
#else
  static MagneticFieldSvc* m_field;
#endif

  /// Fills Q, the field vector
  void fieldGrid( const HepPoint3D& xyz, HepVector3D& fvec ) const;
  void fieldGrid_TE( const HepPoint3D& xyz, HepVector3D& fvec ) const;

  StatusCode parseFile( );       ///< Reads the field map from file
  StatusCode parseFile_TE( );       ///< Reads the field map from file

  std::string path;             // path to data files

  std::string m_filename;        ///< Magnetic field file name
  std::string m_filename_TE;        ///< Magnetic field file name
  //Not open map files for every run
  std::string former_m_filename_TE;//("First Run");
  std::string former_m_filename;//="First Run";

  int m_runmode;                 ///< Run mode
  int m_gridDistance;            ///< grid distance of field map
  int m_outlevel;
  double m_scale;

  double m_Cur_SCQ1_55;
  double m_Cur_SCQ1_89;
  double m_Cur_SCQ2_10;

  std::vector<double> m_Q;       ///< Field vector
  std::vector<double> m_P;       ///< Grid position
  std::vector<double> m_Q_TE;       ///< Field vector
  std::vector<double> m_P_TE;       ///< Grid position
  std::vector<double> m_Q_1;       ///< Field vector
  std::vector<double> m_P_1;       ///< Grid position
  std::vector<double> m_Q_2;       ///< Field vector
  std::vector<double> m_P_2;       ///< Grid position
  double m_Dxyz[3];              ///< Steps in x, y and z
  int    m_Nxyz[3];              ///< Number of steps in x, y and z
  double m_max_FL[3];
  double m_min_FL[3];
  //for tof and emc
  double m_Dxyz_TE[3];              ///< Steps in x, y and z
  int    m_Nxyz_TE[3];              ///< Number of steps in x, y and z
  double m_max_FL_TE[3];
  double m_min_FL_TE[3];

  double m_zOffSet;              ///< The z offset
  double m_zOffSet_TE;              ///< The z offset
  MucMagneticField* m_Mucfield;
  double m_zfield;
  bool m_ifRealField;

  //property
  bool m_useDB;

  bool m_turnOffField;
  bool m_uniField;
public:
  bool m_readOneTime;
  int           m_runFrom;
  int           m_runTo;
  int runNo;

#ifndef BEAN
  IDataProviderSvc* m_eventSvc;
#endif
  //database
  FieldDBUtil::ConnectionDB* m_connect_run;
  std::map<int, std::vector<double> > m_mapMagnetInfo;
  std::map<int, std::vector<double> > m_mapBeamEnergy;
  std::vector<double> beamEnergy;
  std::vector<double> current;
};

#endif  // MAGNETICFIELDSVC_H
