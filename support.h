#ifndef SUPPORT_H
#define SUPPORT_H

#include <stdexcept>
#include <vector>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <iostream>

namespace Sup
{

	double sq(double x);

	const double Pi = 3.14159265358979323846;

	// numerical integration using Simpson's rule
	template <class Fun>
	double integrate(Fun f, double start, double end, unsigned int nr_steps = 10000u)
	{
		if (end < start)
			throw std::runtime_error("integrate: end < start");
		if (end == start)
			return 0.0;

		double I{ 0.0 };
		double stepLength{ (end - start) / double(nr_steps) };

		for (double x{ start }; x <= end; x += stepLength)
			I += (f(x) + f(x + stepLength)) * 0.5 * stepLength;

		return I;
	}

	struct Exp			// parameterised exponent
	{
		double A, k;

		double operator()(double x) const;
		Exp operator * (Exp rhs);
	};

	Exp sq(Exp e);		// exponent squared
	Exp diff(Exp e);	// first order derivative
	Exp diff_2(Exp e);	// second order derivative
	double integrate(Exp e, double start, double end);	// exponent integral exact formula

	struct Pow			// parameterised pow() function
	{
		double power;

		double operator()(double x) const;
	};

	double integrate(Pow p, double start, double end);	// pow() integral exact formula

	struct Line		// linear function
	{
		double b, k;

		double operator()(double x) const;
	};
	
	struct DataPoint	// plot point {x,y}
	{
		double x, y;

		DataPoint() : x{ 0.0 }, y{ 0.0 }
		{ }
		DataPoint(double x_, double y_) : x{ x_ }, y{ y_ }
		{ }
	};

	using DataPoints = std::vector<DataPoint>;	// plot points

	//std::ostream& operator << (std::ostream& str_out, Sup::DataPoints dp)
	//{
	//	for (const auto& p : dp)
	//		str_out << p.x << ' ' << p.y << '\n';
	//}

	std::tuple<double, double> mean_XY(const DataPoints&);				// X and Y mean values
	double cov_XY(const DataPoints&, double mean_x, double mean_y);		// two samples' covariance
	double var_X(const DataPoints&, double mean_x);						// variance 
	Line RegressionLine(const DataPoints&);								// linear regression

	DataPoints logY(DataPoints dp);										// log of a vector
	Exp fitExp(DataPoints);												// least squares exponent fitting


	struct Spline			// "stitching" function
	{
		double a, b;		// domain
		double k_a, k_b;	// special params

		Spline(double a, double b, double k_a, double k_b);
		Line l;
		Exp e;
		double operator()(double x) const;
	};

	template <class Fun>
	DataPoints makeSamples(Fun f, double start, double end, unsigned samples_nr = 1000u)
	{
		double step{ (end - start) / samples_nr };
		DataPoints pts(samples_nr);
		int i{ 0 };

		std::generate_n(begin(pts), pts.size(), [&](){ 
			auto x = start + step * i;
			auto y = f(x);
			i++;
			return DataPoint(x, y);
		});

		return pts;
	}

	template<class T>
	typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
		close(T x, T y, int ulp)
	{
		// the machine epsilon has to be scaled to the magnitude of the values used
		// and multiplied by the desired precision in ULPs (units in the last place)
		return std::abs(x - y) < std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
			// unless the result is subnormal
			|| std::abs(x - y) < std::numeric_limits<T>::min();
	}

} // Sup namespace


#endif // SUPPORT_H