#include <iomanip>
#include <omp.h>

#include "precision.h"
#include "linsol_error.h"

linsol_error::linsol_error(std::string const & __restrict__ function_name,
                           std::string const & __restrict__ error_message) {

   std::cout << "FUNCTION: " << function_name << " - ERROR: " << error_message << std::endl;

}




linsol_error::linsol_error(std::string const & __restrict__ function_name,
                           std::string const & __restrict__ error_message,
                           unsigned long int nbytes_requested) {

   #pragma omp critical (error_print)
   {


      type_OMP_iReg mythid = omp_get_thread_num();


      std::cout << "THREAD ID: " << mythid << std::endl;
      std::cout << "FUNCTION: " << function_name << " - ERROR: "<< error_message << std::endl;
      rExt req = static_cast<rExt>(nbytes_requested) / static_cast<rExt>(MemUM);
      std::cout << std::fixed << std::setprecision(5);
      std::cout << "REQUESTED:           " << req << " " << UM_label << std::endl;
      std::cout.flush();

   }

}
