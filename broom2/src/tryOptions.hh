//---------------------------------------------------------------------------
#ifndef __TRYOPTIONS_HH__
#define __TRYOPTIONS_HH__

#include <string>
#include <vector>
#include <cstdlib>

#include "mfem.hpp"

/// Look for a particular option on the command line, and set the value.  Pass
/// in the default int.
inline int tryOption(const std::string option,
                     const std::string& message,
                     int argc,
                     char* const argv[],
                     const int defaultVal)
{
   int value = defaultVal;
   for (int i = 1; i < argc; ++i)
   {
      if (option == argv[i])
      {
         if (i + 1 < argc)
         {
            char* end;
            int temp = strtol(argv[i + 1], &end, 10);
            if (end != argv[i + 1])
            {
               value = temp;
               std::cout << option << " = " << value << "\n        " << message
                         << " (default is " << defaultVal << ")\n";
               return value;
            }
            else if (argv[i + 1][0] == '-')
            {
               MFEM_ABORT(option << " option requires a int parameter");
            }
            else
            {
               MFEM_ABORT("Unable to parse option "
                          << option << "'s required int parameter.  You specified "
                          << argv[i + 1]);
            }
         }
         else
         {
            MFEM_ABORT(option << " option requires a int parameter");
         }
      }
   }
   std::cout << option << " = " << value << "\n        " << message
             << " (using default)\n";
   return value;
}

//---------------------------------------------------------------------------

/// Look for a particular option on the command line, and set the value.  Pass
/// in the default double.
inline double tryOption(const std::string option,
                        const std::string& message,
                        int argc,
                        char* const argv[],
                        const double defaultVal)
{
   double value = defaultVal;
   for (int i = 1; i < argc; ++i)
   {
      if (option == argv[i])
      {
         if (i + 1 < argc)
         {
            char* end;
            double temp = strtod(argv[i + 1], &end);
            if (end != argv[i + 1])
            {
               value = temp;
               std::cout << option << " = " << value << "\n        " << message
                         << " (default is " << defaultVal << ")\n";
               return value;
            }
            else if (argv[i + 1][0] == '-')
            {
               MFEM_ABORT(option << " option requires a double parameter");
            }
            else
            {
               MFEM_ABORT("Unable to parse option "
                          << option << "'s required double parameter.  You specified "
                          << argv[i + 1]);
            }
         }
         else
         {
            MFEM_ABORT(option << " option requires a double parameter");
         }
      }
   }
   std::cout << option << " = " << value << "\n        " << message
             << " (using default)\n";
   return value;
}

//---------------------------------------------------------------------------

/// Look for a particular option on the command line, and set the value.  Pass
/// in the default string.
inline std::string tryOption(const std::string option,
                             const std::string& message,
                             int argc,
                             char* const argv[],
                             const std::string& defaultVal)
{
   std::string value = defaultVal;
   for (int i = 1; i < argc; ++i)
   {
      if (option == argv[i])
      {
         if (i + 1 < argc)
         {
            if (argv[i + 1][0] == '-')
            {
               MFEM_ABORT(option << " option requires a string parameter");
            }
            else
            {
               value = argv[i + 1];
               std::cout << option << " = " << value << "\n        " << message
                         << " (default is " << defaultVal << ")\n";
               return value;
            }
         }
         else
         {
            MFEM_ABORT(option << " option requires a string parameter");
         }
      }
   }
   std::cout << option << " = " << value << "\n        " << message
             << " (using default)\n";
   return value;
}

//---------------------------------------------------------------------------

/// Look for a particular option on the command line, and set the value. Changes
/// the boolean to true.
inline bool tryOption(const std::string trueOption,
                      const std::string falseOption,
                      const std::string& message,
                      int argc,
                      char* const argv[],
                      const bool defaultVal)
{
   bool value = defaultVal;
   for (int i = 1; i < argc; ++i)
   {
      if (trueOption == argv[i])
      {
         value = true;
         std::cout << trueOption << " or " << falseOption << " = " << trueOption
                   << "\n        " << message << " (default is "
                   << (defaultVal ? trueOption : falseOption) << ")\n";
         return value;
      }
      if (falseOption == argv[i])
      {
         value = false;
         std::cout << trueOption << " or " << falseOption << " = " << falseOption
                   << "\n        " << message << " (default is "
                   << (defaultVal ? trueOption : falseOption) << ")\n";
         return value;
      }
   }
   std::cout << trueOption << " or " << falseOption << " = "
             << (defaultVal ? trueOption : falseOption) << "\n        " << message
             << " (using default)\n";
   return value;
}

#endif  // __TRYOPTIONS_HH__
