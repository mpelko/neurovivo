#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," NEURON_stuff/HayStuff/mod/CaDynamics_E2.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/Ca_HVA.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/Ca_LVAst.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/Ih.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/Im.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/K_Pst.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/K_Tst.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/NaTa_t.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/NaTs2_t.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/Nap_Et2.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/SK_E2.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/SKv3_1.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/epsp.mod");
    fprintf(stderr," NEURON_stuff/HayStuff/mod/vecevent.mod");
    fprintf(stderr, "\n");
  }
  _CaDynamics_E2_reg();
  _Ca_HVA_reg();
  _Ca_LVAst_reg();
  _Ih_reg();
  _Im_reg();
  _K_Pst_reg();
  _K_Tst_reg();
  _NaTa_t_reg();
  _NaTs2_t_reg();
  _Nap_Et2_reg();
  _SK_E2_reg();
  _SKv3_1_reg();
  _epsp_reg();
  _vecevent_reg();
}
