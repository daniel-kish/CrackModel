#ifndef CRACK_GROWTH_MODEL_H
#define CRACK_GROWTH_MODEL_H

#include "support.h"

class CrackGrowthModel
{
public:
	CrackGrowthModel(double l0 = 0.005, double f = 1.0);


	void nextState();
	void simulate();
	double lenMaxPermissible();
	double time() const { return t; }

	double l = 0.005;	// текущая длина трещины

	double l0 = 0.005;	// начальная длина трещины
	double t = 0.0;
	bool isEmbrittlement = false;

	double Dt_e, Dt_f;

	double a;
	double b;

	double C0;
	double mean_Ca;

	double C = 1.095E-12;
	double n = 3.24;
	double f = 0.0;

	double Omega = 2.0;
	double D = 1.0E-10;
	double V_H = 1.96E-6;
	double R = 8.314;
	double T = 293.0;

	Sup::Exp phi;		// текущая координатная функция (м-д Галеркина)

	double delta = 0.0002;
	double a0 = 1.0E-5;
	
	double B = 10.0;
	double K_max = 80.0;
	double K_min = 10.0;
	double alpha_1 = 2.0, beta_1 = 2.0;
	double sigma_out = 140.0;
	double alpha_2 = 2.0, beta_2 = 2.0;
	double Y = 1.0; // К-тарировка

	double k_a = 0.0025, k_b = 0.000001;

	double preruptureAreaLen();
	
	double lenFracture();
	double timeToRupture_embrittlement();
	double timeToRupture_fatigue();
	double K_I();
	double criticalConc();
	double coef();
};

#endif // CRACK_GROWTH_MODEL_H