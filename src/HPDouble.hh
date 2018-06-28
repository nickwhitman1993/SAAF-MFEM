#ifndef __HPDOUBLE_HH__
#define __HPDOUBLE_HH__

/// \file HPDouble.hh
/// \brief Do a sum of many doubles with extra precision.

/*! \brief A high precision double, constructed from the double-double library
 * of Hida et al.
 *
 *  This class will sum numbers with extra precision.  It isn't meant to
 *  be a double-double library, so any time you want to do a multiply,
 *  divide, or some special function, you need to convert back to a double.
 *
 *  The algorithm uses two doubles, the main sum and a correction.  Sums are
 *  evaluated in steps using four different tools derived directly from the
 *  doubl-double library of Hida et al. which inturn were derived from the work
 *  of Shawchuk on arbitary precission arithmatic.
 *
 * double quick_two_sum(a,b,c)
 *  \f{align}{
 *      s &= a+b\\
 *      c &= b-(s-a)\\
 *      return s
 *  \f}
 *
 * double two_sum(a,b,c)
 *  \f{align}{
 *      s &= a+b\\
 *      c &= b-(s-a)\\
 *      return s
 *  \f}
 *
 * Adding a double-double (a0,a1) value to a double (b0) can be done by:
 *  \f{align}{
 *      temp_s &= two_sum(a0,b0,temp_c)\\
 *      temp_c +&= a1\\
 *      s &= quick_two_sum(temp_s,temp_c,c)
 *      return s
 *  \f}
 * Similarly adding a double-double (a0,a1) value to a double-double (b0,b1) can
 * be done by:
 *  \f{align}{
 *      temp_s &= two_sum(a0,b0,temp_c)\\
 *      temp_c +&= a1\\
 *      temp_c +&= b1\\
 *      s &= quick_two_sum(temp_s,temp_c,c)
 *      return s
 *  \f}
 *
 *  There might be trouble in an x86 environment, in it isn't exactly following
 *  64-bit IEEE-754 rules.  (The registers are 80 bit.)  So far, we haven't seen
 *  anything.  It may be possible to set some CPU control flags, but we don't.
 *  See the qd library http://crd.lbl.gov/~dhbailey/mpdist/ for examples.

 *  References:
 *
 *  - Y. Hida, X. S. Li, and D. H. Bailey, "Library for Double-Double and
 *  Quad-Double Arithmetic," December 29, 2007, (Berkley National Lab report)
 *
 *  - J. R. Shewchuk, "Adaptive Precision Floating-Point Arithmetic and Fast
 *  Robust Geometric Predicates," Discrete & Computational Geometry
 *  18(3):305-363, October 1997.
 *
 */

class HPDouble
{
  public:
   //---------------------------------------------------

   /// Default constructor sets everything to zero.
   HPDouble() : mValue(0.0), mCorrection(0.0) {}
   //---------------------------------------------------

   /// Allow implicit conversion from a double
   HPDouble(const double value) : mValue(value), mCorrection(0.0) {}
   //---------------------------------------------------

   /// Set both parameters.  Useful for implementing operator-().
   HPDouble(const double value, const double correction)
       : mValue(value), mCorrection(correction)
   {
   }

   //---------------------------------------------------

   /// Copy constructor
   HPDouble(const HPDouble& rhs) : mValue(rhs.mValue), mCorrection(rhs.mCorrection) {}
   //---------------------------------------------------

   /// Allow implicit conversion to a double:  Is this bad?
   operator double() const
   {
      double temp = mCorrection;
      temp += mValue;
      return temp;
   }

   //---------------------------------------------------

   /// Set from a double.
   const HPDouble& operator=(const double& x)
   {
      mValue = x;
      mCorrection = 0.0;
      return *this;
   }

   //---------------------------------------------------

   /// Set from an HPDouble
   const HPDouble& operator=(const HPDouble& rhs)
   {
      mValue = rhs.mValue;
      mCorrection = rhs.mCorrection;
      return *this;
   }

   //---------------------------------------------------

   /// Expansion_sum as defined by Shewchuk that only carries two doubles
   // The lower order corrections are added to the highest correction
   HPDouble& operator+=(const HPDouble& rhs)
   {
      double s, e;
      s = two_sum(mValue, rhs.mValue, e);
      e += mCorrection;
      e += rhs.mCorrection;
      mValue = quick_two_sum(s, e, mCorrection);
      return *this;
   }

   //---------------------------------------------------

   /// The same as above, but no error term to add in.
   HPDouble& operator+=(const double rhs)
   {
      double s1, s2;
      s1 = two_sum(mValue, rhs, s2);
      s2 += mCorrection;
      mValue = quick_two_sum(s1, s2, mCorrection);
      return *this;
   }

   //---------------------------------------------------

   /// Subtraction is just addition with the right minus signs
   HPDouble& operator-=(const HPDouble& rhs)
   {
      double s, e;
      s = two_diff(mValue, rhs.mValue, e);
      e += mCorrection;
      e -= rhs.mCorrection;
      mValue = quick_two_sum(s, e, mCorrection);
      return *this;
   }

   //---------------------------------------------------

   /// Same as subtracting an HP double, just with a minus sign or two.
   HPDouble& operator-=(const double rhs)
   {
      double s1, s2;
      s1 = two_diff(mValue, rhs, s2);
      s2 += mCorrection;
      mValue = quick_two_sum(s1, s2, mCorrection);
      return *this;
   }

   //---------------------------------------------------

   /// Unary + operator
   const HPDouble& operator+() const { return *this; }
   //---------------------------------------------------

   /// Unary - operator
   const HPDouble operator-() const
   {
      HPDouble tmp(-mValue, -mCorrection);
      return tmp;
   }

   //---------------------------------------------------

   /// Access to the value directly, handy for ostream operator, if we had one.
   double getValue() const { return mValue; }
   //---------------------------------------------------

   /// Access to the correction directly, handy for ostream operator, if we had
   /// one.
   double getCorrection() const { return mCorrection; }
   //---------------------------------------------------
   //---------------------------------------------------

  private:
   /// The current value of the double
   double mValue;
   /// The correction that's used to improve the next summation.
   double mCorrection;
   /// Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|
   inline double quick_two_sum(double a, double b, double& err)
   {
      double s = a + b;
      err = b - (s - a);
      return s;
   }

   /// Computes fl(a-b) and err(a-b).  Assumes |a| >= |b|
   inline double quick_two_diff(double a, double b, double& err)
   {
      double s = a - b;
      err = (a - s) - b;
      return s;
   }

   /// Computes fl(a+b) and err(a+b)
   inline double two_sum(double a, double b, double& err)
   {
      double s = a + b;
      double bb = s - a;
      err = (a - (s - bb)) + (b - bb);
      return s;
   }

   /// Computes fl(a-b) and err(a-b)
   inline double two_diff(double a, double b, double& err)
   {
      double s = a - b;
      double bb = s - a;
      err = (a - (s - bb)) - (b + bb);
      return s;
   }
};

//----------------------------------------------------------------
//----------------------------------------------------------------

/// \name Binary Operators for HPDouble
/// We are only defining the addition and subtraction operators
/// for HPDouble.  All the other operators rely on the implicit conversion
/// to a double.  Doing more essentially makes this a quad-precision
/// class, and we want to avoid that.
///@{

inline const HPDouble operator+(const HPDouble& lhs, const HPDouble& rhs)
{
   // Enable named return value optimization
   HPDouble tmp(lhs);
   tmp += rhs;
   return tmp;
}

//----------------------------------------------------------------

inline const HPDouble operator+(const HPDouble& lhs, const double rhs)
{
   // Enable named return value optimization
   HPDouble tmp(lhs);
   tmp += rhs;
   return tmp;
}

//----------------------------------------------------------------

inline const HPDouble operator+(const double lhs, const HPDouble& rhs)
{
   // Enable named return value optimization
   HPDouble tmp(rhs);
   tmp += lhs;
   return tmp;
}

//----------------------------------------------------------------

inline const HPDouble operator-(const HPDouble& lhs, const HPDouble& rhs)
{
   // Enable named return value optimization
   HPDouble tmp(lhs);
   tmp -= rhs;
   return tmp;
}

//----------------------------------------------------------------

inline const HPDouble operator-(const HPDouble& lhs, const double rhs)
{
   // Enable named return value optimization
   HPDouble tmp(lhs);
   tmp -= rhs;
   return tmp;
}

//----------------------------------------------------------------

inline const HPDouble operator-(const double lhs, const HPDouble& rhs)
{
   // Enable named return value optimization
   HPDouble tmp(-rhs);
   tmp += lhs;
   return tmp;
}
///@}

//----------------------------------------------------------------

#endif  // __HPDOUBLE_HH__
