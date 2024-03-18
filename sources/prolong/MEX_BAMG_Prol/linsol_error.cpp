#include <iomanip>
#include <omp.h>
using namespace std;

#include "precision.h"
#include "linsol_error.h"




linsol_error::linsol_error(const string __restrict__ function_name,
                           const string __restrict__ error_message) {

   cout << "FUNCTION: " << function_name << " - ERROR: " << error_message << endl;

}




linsol_error::linsol_error(const string __restrict__ function_name,
                           const string __restrict__ error_message,
                           unsigned long int nbytes_requested) {

   #pragma omp critical (error_print)
   {


      type_OMP_iReg mythid = omp_get_thread_num();


      cout << "THREAD ID: " << mythid << endl;
      cout << "FUNCTION: " << function_name << " - ERROR: "<< error_message << endl;
      rExt req = static_cast<rExt>(nbytes_requested) / static_cast<rExt>(MemUM);
      cout << fixed << setprecision(5);
      cout << "REQUESTED:           " << req << " " << UM_label << endl;
      cout.flush();

   }

}
