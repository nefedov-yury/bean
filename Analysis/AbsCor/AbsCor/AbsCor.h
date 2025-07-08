#ifndef Analysis_AbsCor_H
#define Analysis_AbsCor_H

#include <string>
#include "ReadDst.h"
#include "TGraph2DErrors.h"

//--------------------------------------------------------------------
// This class is a replacement for class AbsCor from the Boss:
// 'boss-version'/Analysis/PhotonCor/AbsCor
//--------------------------------------------------------------------

class AbsCor {

   public:
      AbsCor(const std::string& path,
            bool UseTof = true,
            bool DoDataCor = true
      );

      // read corrections from given files:
      void ReadDatac3p(std::string DataPathc3p, bool info=true);
      void ReadParMcCor(std::string paraPath, bool info=true);
      void ReadCorFunpara(std::string CorFunparaPath, bool info=true);

      // apply corrections:
      void AbsorptionCorrection(ReadDst* selector);
      void SuppressHotCrystals(ReadDst* selector);

      // functions to managing parameters of algorithm:
      void UseTof(bool use = true) {
         usetof = use;
      }
      void DoDataCor(bool cor = true) {
         dodatacor = cor;
      }
      void SuppressHotCrystals(bool use = false) {
         hotcellmask = use;
      }
      void DoPi0Cor(bool cor = true) {
         dopi0Cor = cor;
      }
      void McUseTof(bool use = true) {
         MCuseTof = use;
      }
      void SetMCCorUseFunction(bool use = true) {
         MCCorUseFunction = use;
      }

   private:
      // parameters of AbsCor:
      bool          usetof;
      bool          dodatacor;

      bool          hotcellmask;
      bool          dopi0Cor;
      bool          MCuseTof;

      bool          MCCorUseFunction;
      std::string   m_DataPathc3ptof;
      std::string   m_CorFunparaPath;

      // arrays and structures for internal needs:
      double        ai[4];

      int           hrunstart[10];
      int           hrunend[10];
      int           hcell[10];

      double        e25min[28];
      double        e25max[28];

      double        m_corFunPar[28][6];

      TGraph2DErrors* dt;       // Shower energy correction
      TGraph2DErrors* dtErr;    // Energy error

      // correction functions:
      double ECorrMC(double eg, double theid) const;
      double ErrMC(double eg, double theid) const;
      double E25min(int n) const {
         return e25min[n];
      }
      double E25max(int n) const {
         return e25max[n];
      }

      double ECorrFunctionMC(double eg, double theid) const;
};
#endif
