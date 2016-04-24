#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>

#include "support.h"
#include "CrackGrowthModel.h"

double time_fat(double l, double l0, double f)
{
	using Sup::Pi;
	const double n = 3.24, C = 1.095E-12, sigma = 140.0;

	double nom = 2.0 * ( pow(l0, 1-n/2) - pow(l, 1-n/2) );
	double denom = pow(Pi,n/2)*pow(sigma,n)*C*(n-2);
	return nom / denom / f;
}


double get_l0_cr(double f, double l0min)
{
	using namespace Sup;

	auto life_emb = [&](double l0)	{
		CrackGrowthModel m(l0, f);

		try	{ m.simulate(); }

		catch (const std::runtime_error& e)	{
			std::cerr << e.what() << '\n';
		}

		return m.time();
	};
	CrackGrowthModel m;
	double Lfr = m.lenMaxPermissible();

	auto life_fat = [&](double l0){ return time_fat(Lfr, l0, f); };

	auto close = [](double x, double y, double eps)	{
		return abs(x - y) < eps;
	};
	
	double step = (Lfr - l0min) / 300.0;
	double max{ life_fat(l0min) }, min{ life_emb(l0min) }, eps{ 250.0 };

	for (double l0 = l0min; l0 <= Lfr; l0 += step)
	{
		if (close(life_emb(l0), life_fat(l0), eps))
			return l0;
	}
}

std::string fileNameNo(int nr)
{
	std::ostringstream name;
	name << "file_" << nr;
	return name.str();
}

using namespace std;
int main(int argc, char* argv[])
{
	CrackGrowthModel m;
	
	double L_max = m.lenMaxPermissible();

	ofstream outfile("myfile.dat");
	while (m.l <= L_max)
	{
		std::cout << m.l << ' ' << m.t << ' ' << m.phi.A << ' ' << m.phi.k << '\n';
		m.nextState();
	}
}
