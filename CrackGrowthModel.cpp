#include "CrackGrowthModel.h"
#include <cmath>
#include "support.h"
#include <iostream>
#include <fstream>

CrackGrowthModel::CrackGrowthModel(double l0_, double f_) : l0{ l0_ }, f{ f_ }
{
	l = l0;
	t = 0.0;
	k_a = 0.0025;
	k_b = 0.000001;

	a = a0;
	b = 500 * a;

	phi = Sup::fitExp(Sup::makeSamples(Sup::Spline(a, b, k_a, k_b), 0.0, b, 2500u));

	mean_Ca = Sup::integrate(phi, 0.0, a) / a;

	Dt_e = timeToRupture_embrittlement();
	Dt_f = timeToRupture_fatigue();

	if (Dt_e < Dt_f) {
		C0 = std::pow((K_I() - K_min) / (K_max - K_min), alpha_1);
		C0 = std::pow(1.0 - C0, -1.0 / beta_1) * Omega * mean_Ca;
	}
	else 
		C0 = exp(-coef() * Dt_f);
	auto m = 1.0 / C0;
}

void CrackGrowthModel::nextState()
{ 
	auto a_prev = a;
	l += a;
	a = preruptureAreaLen();
	b = 500.0*a;

	k_a = phi(a_prev + a) / C0;
	k_b = phi(a_prev + b) / C0;

	if (Dt_e < Dt_f) {
		isEmbrittlement = true;
		t += Dt_e;
	}
	else {
		isEmbrittlement = false;
		t += Dt_f;
	}

	phi = Sup::fitExp(Sup::makeSamples(Sup::Spline(a, b, k_a, k_b), 0.0, b, 2500u));

	mean_Ca = Sup::integrate(phi, 0.0, a) / a;

	Dt_e = timeToRupture_embrittlement();
	Dt_f = timeToRupture_fatigue();

	if (Dt_e < Dt_f) {
		C0 = std::pow((K_I() - K_min) / (K_max - K_min), alpha_1);
		C0 = std::pow(1.0 - C0, -1.0 / beta_1) * Omega * mean_Ca;
	}
	else
		C0 = exp(-coef() * Dt_f);
}

void CrackGrowthModel::simulate()
{
	std::ofstream out("C:\\Users\\Daniel\\Documents\\Visual Studio 2013\\Projects\\StochasticBatch\\Release\\l0cr_f\\newfolder\\new.dat");
	int N{ 0 };
	auto L_max = lenMaxPermissible();
	while (l <= L_max)
	{
		out << l << ' ' << t << ' ' << mean_Ca << '\n';
		nextState();
		N++;
	}
}

double CrackGrowthModel::preruptureAreaLen()
{
	double L_max = lenMaxPermissible();
	double A = std::pow((L_max - l) / (L_max - l0), beta_2);
	A = std::pow(1.0 - A, 1.0 / alpha_2);
	return a0*(A*(B - 1.0) + 1.0);
}

double CrackGrowthModel::lenMaxPermissible()
{
	return (1.0 - delta) * lenFracture();
}

double CrackGrowthModel::lenFracture()
{
	using Sup::sq;
	using Sup::Pi;

	return sq(K_max / (sigma_out * Y)) / Pi;
}

double CrackGrowthModel::timeToRupture_embrittlement()
{
	using Sup::integrate;
	using Sup::sq;
	using Sup::diff_2;
	using Sup::Pi;
	using Sup::Exp;

	double a = preruptureAreaLen();

	double Beta = (K_I() - K_min) / (K_max - K_min);
	Beta = 1.0 - std::pow(Beta, alpha_1);
	Beta = std::pow(Beta, 1.0 / beta_1);

	Beta *= a / Omega;
	Beta /= integrate(phi, 0.0, a);

	double dzeta = (D*V_H) / (3.0*R*T);

	double num = integrate(sq(phi), a, b);
	double denom = D * integrate(diff_2(phi)*phi, a, b);

	Exp g = diff(phi) * phi;
	auto integrand = [&] (double x)	{	return g(x) / (x*sqrt(x));	};

	double integ = integrate(integrand, a, b, 10000u) / (2.0*sqrt(Pi)); 
	double K_IVal = K_I();

	denom -= dzeta * 1.0e6 /* ??? */  * K_IVal * integ;
	
	double Gamma = num / denom;
	double DeltaT = Gamma * std::log(Beta);
	return DeltaT;
}

double CrackGrowthModel::timeToRupture_fatigue()
{
	using Sup::Pi;
	using std::pow;
	using Sup::integrate;

	double a = preruptureAreaLen();

	double power = -n*0.5;
	Sup::Pow p{ power };

	double integ = integrate(p, l, l + a);
	integ /= (C * pow(sigma_out,n) * pow(Pi,n*0.5));
	
	return integ / f;
}

double CrackGrowthModel::K_I()
{
	using Sup::Pi;
	return sigma_out * std::sqrt(Pi*l) * Y;
}

double CrackGrowthModel::criticalConc()
{
	double Beta = (K_I() - K_min) / (K_max - K_min);
	Beta = 1.0 - std::pow(Beta, alpha_1);
	Beta = std::pow(Beta, 1.0 / beta_1);
	Beta *= C0/Omega;

	return Beta;
}

double CrackGrowthModel::coef()
{
	using namespace Sup;

	double dzeta = (D*V_H) / (3.0*R*T);

	double denom = integrate(sq(phi), a, b);
	double num = D * integrate(diff_2(phi)*phi, a, b);

	Exp g = diff(phi) * phi;
	auto integrand = [&](double x)	{	return g(x) / (x*sqrt(x));	};

	double integ = integrate(integrand, a, b, 10000u) / (2.0*sqrt(Pi));
	double K_IVal = K_I();

	num -= dzeta * 1.0e6 * K_IVal * integ;

	return num / denom;
}