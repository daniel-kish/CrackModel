#include "support.h"
#include <tuple>

double Sup::sq(double x)
{
	return x*x;
}

double Sup::Exp::operator()(double x) const
{
	return A * std::exp(k * x);
}


Sup::Exp Sup::sq(Sup::Exp e)
{
	Sup::Exp e2{ Sup::sq(e.A), 2.0*e.k };
	return e2;
}

Sup::Exp Sup::diff(Sup::Exp e)
{
	Sup::Exp e2{ e.A*e.k, e.k };
	return e2;
}

Sup::Exp Sup::diff_2(Sup::Exp e)
{
	return diff(diff(e));
}

Sup::Exp Sup::Exp::operator* (Sup::Exp rhs)
{
	Sup::Exp res{ rhs.A*A, rhs.k + k };
	return res;
}

double Sup::integrate(Sup::Exp e, double start, double end)
{
	using std::exp;

	if (end < start)
		throw std::runtime_error("integrate: end < start");

	double integralValue = exp(start*e.k) - exp(end*e.k);
	integralValue *= -e.A / e.k;
	return integralValue;
}

double Sup::Pow::operator()(double x) const
{
	return std::pow(x, power);
}


double Sup::integrate(Sup::Pow p, double start, double end)
{
	double I = 1.0 / (p.power + 1.0);
	I *= (std::pow(end, p.power + 1.0) - std::pow(start, p.power + 1.0));
	return I;
}

double Sup::Line::operator()(double x) const
{
	return b + k*x;
}

std::tuple<double, double> Sup::mean_XY(const Sup::DataPoints& dp)
{
	double mean_x{}, mean_y{};

	for (const auto& p : dp)
	{
		mean_x += p.x;
		mean_y += p.y;
	}
	mean_x /= dp.size();
	mean_y /= dp.size();

	return std::make_tuple(mean_x, mean_y);
}

double Sup::cov_XY(const Sup::DataPoints& dp, double mean_x, double mean_y)
{
	double cov{};

	for (const auto& p : dp)
		cov += (p.x - mean_x) * (p.y - mean_y);

	return cov / dp.size();
}

double Sup::var_X(const Sup::DataPoints& dp, double mean_x)
{
	double var{};

	for (const auto& p : dp)
		var += sq( (p.x - mean_x) );

	return var / dp.size();
}

Sup::Line Sup::RegressionLine(const Sup::DataPoints& dp)
{
	double mean_x{}, mean_y{};

	std::tie(mean_x, mean_y) = mean_XY(dp);

	double cov = cov_XY(dp, mean_x, mean_y);
	double var = var_X(dp, mean_x);

	double slope{ cov / var };

	return Sup::Line{ mean_y - slope*mean_x, slope};
}

Sup::DataPoints Sup::logY(Sup::DataPoints dp)
{
	// логарифмируем Y каждой точки
	for (auto& point : dp)
		point.y = std::log(point.y);

	return dp;
}

Sup::Exp Sup::fitExp(Sup::DataPoints dp)
{
	Line l = Sup::RegressionLine(Sup::logY(move(dp)));

	return Sup::Exp{ std::exp(l.b), l.k };
}

Sup::Spline::Spline(double a, double b, double k_a, double k_b)
	: a{ a }, b{ b }, k_a{ k_a }, k_b{ k_b }
{	}

double Sup::Spline::operator()(double x) const
{
	if (x < 0.0 || x > b)
		throw std::runtime_error("spline arguments out of [a; b]");

	Line l{ 1.0, -(1.0 - k_a)/a };
	Exp e{ k_a, std::log(k_b / k_a) };

	if ( close(x,a,16) )
		return k_a;
	else if ( close(x,a,16) )
		return b;

	if (x >= 0.0 && x <= a)
		return l(x);
	else
		return e((x - a) / (b - a));
}