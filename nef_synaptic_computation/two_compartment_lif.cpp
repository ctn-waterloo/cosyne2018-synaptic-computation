/*
 *   Conductance based synapses in Nengo
 *   Copyright (C) 2018  Andreas Stöckel
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>

#pragma pack(push, 8)
struct Params {
	double gC = 100e-9;
	double gLDen = 12.5e-9;
	double gLSom = 50e-9;
	double cmDen = 0.25e-9;
	double cmSom = 1.0e-9;
	double EE = 20e-3;
	double EI = -75e-3;
	double EL = -65e-3;
	double v_thresh = -50e-3;
	double tau_ref = 2e-3;
};

struct alignas(32) fmat_2x2 {
	double a, b, c, d;
};
#pragma pack(pop)

static inline fmat_2x2 inv(const fmat_2x2 &in)
{
	const double &a = in.a, &b = in.b, &c = in.c, &d = in.d;
	const double det = a * d - b * c;

	return fmat_2x2{d / det, -b / det, -c / det, a / det};
}

static inline fmat_2x2 expm(const fmat_2x2 &in, double t)
{
	const double &a = in.a, &b = in.b, &c = in.c, &d = in.d;
	const double D = std::sqrt(a * a + d * d - 2.0f * a * d + 4.0f * b * c);

	const double lambda_1 = 0.5f * (d + a + D);
	const double lambda_2 = 0.5f * (d + a - D);

	const double c_1 = (2.0f * c / (d - a + D));
	const double c_2 = (2.0f * c / (d - a - D));
	const double c_d = c_1 - c_2;

	const double e_1 = std::exp(lambda_1 * t);
	const double e_2 = std::exp(lambda_2 * t);

	return fmat_2x2{(c_1 * e_1 - c_2 * e_2) / c_d, (e_1 - e_2) / c_d,
	                (-c_1 * e_1 * c_2 + c_2 * e_2 * c_1) / c_d,
	                (-c_2 * e_1 + c_1 * e_2) / c_d};
}

extern "C" {
void make_noise(double tau, double rate, double dt, double *xs, uint64_t n)
{
	std::default_random_engine gen;
	std::uniform_real_distribution<double> d_weights(0.0, 2.0);
	std::exponential_distribution<double> d_spikes(rate);

	const double f1 = 1.0 / (tau * rate);
	const double f0 = (1.0 - dt / tau);
	double t = d_spikes(gen);
	xs[0] = 1.0;
	for (size_t i = 1; i < n; i++) {
		xs[i] = xs[i - 1] * f0;
		t = t - dt;
		while (t < 0) {
			xs[i] += f1 * d_weights(gen);
			t += d_spikes(gen);
		}
	}
}

void simulate_two_compartment_lif_conductances(const Params *params,
                                               const double *gEs,
                                               const double *gIs,
                                               double *state,
                                               double *spikes,
                                               double *i_syn,
                                               uint64_t n, double dt)
{
	const Params p = *params;  // For more convenient access

	// Pre-compute some constants
	const double Ab = p.gC / p.cmDen;
	const double Ac = p.gC / p.cmSom;
	const double Ad = -(p.gLSom + p.gC) / p.cmSom;
	const double b2 = p.EL * p.gLSom / p.cmSom;

	// State variables (v1, v2, tRef)
	double v1 = state ? state[0] : p.EL;
	double v2 = state ? state[1] : p.EL;
	double tRef = state ? state[2] : 0.0;

	for (size_t i = 0; i < n; i++) {
		// Fetch the current input conductances
		const double gE = gEs[i];
		const double gI = gIs[i];

		// Assemble and invert the A matrix
		const fmat_2x2 A{-(gE + gI + p.gLDen + p.gC) / p.cmDen, Ab, Ac, Ad};
		const fmat_2x2 AInv = inv(A);
		const fmat_2x2 AExp = expm(A, dt);

		// Assemble the b vector and calculate the equilibrium potentials
		const double b1 = (p.EE * gE + p.EI * gI + p.EL * p.gLDen) / p.cmDen;
		const double vEq1 = -AInv.a * b1 - AInv.b * b2;
		const double vEq2 = -AInv.c * b1 - AInv.d * b2;

		// Compute the state in time dt
		double v1Next = vEq1 + AExp.a * (v1 - vEq1) + AExp.b * (v2 - vEq2);
		double v2Next = vEq2 + AExp.c * (v1 - vEq1) + AExp.d * (v2 - vEq2);

		// Spike mechanism
		if (v2Next > p.v_thresh) {
			tRef = p.tau_ref + dt;
			spikes[i] = 1.0f / dt;
		}
		if (tRef > dt) {
			tRef = tRef - dt;
			v2Next = p.EL;
			if (i_syn) {
				i_syn[i] = std::numeric_limits<double>::quiet_NaN();
			}
		} else if (i_syn) {
			i_syn[i] = (v1Next - v2Next) * p.gC;
		}

		// Assign the next state to the current state
		v1 = v1Next, v2 = v2Next;
	}

	// Write the current state back
	if (state) {
		state[0] = v1;
		state[1] = v2;
		state[2] = tRef;
	}
}
}

