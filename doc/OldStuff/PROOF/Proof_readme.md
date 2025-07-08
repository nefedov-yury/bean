### ROOT PROOF

* The PROOF framework is marked as deprecated in the ROOT since
  version 6.26, and was completely removed in version 6.32.
  Therefore, the BEAN was rewritten so that it can be compiled with
  or without the PROOF support.

* You should use the cmake option `-DUSE_PROOF=ON` to enable the PROOF
  functionality in the BEAN program.
  The default setting is `USE_PROOF=OFF`.

* `make proofbean` command creates PAR files for ROOT-PROOF.

*  The work of BEAN on the PROOF cluster has not been tested for a very
   long time.

