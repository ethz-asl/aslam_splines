/*
 * NumericIntegrator.hpp
 *
 *  Created on: Oct 7, 2012
 *      Author: hannes
 */

#ifndef NUMERICINTEGRATOR_HPP_
#define NUMERICINTEGRATOR_HPP_

#include <cmath>

namespace numeric_integrator {
	template <typename TValue, typename TArgScalar>
	struct Integrand {
		typedef TValue ValueT;
		typedef TArgScalar ArgScalarT;

		typedef ValueT (&IntegrandFunctionT)(TArgScalar & t) ;


		Integrand(IntegrandFunctionT integrand) : _integrand(integrand){}

		inline TValue operator () (const TArgScalar arg) const {
			return _integrand(arg);
		}

	private:
		IntegrandFunctionT _integrand;
	};

	template <typename TValue, typename TArgScalar>
	inline Integrand<TValue, TArgScalar> createIntegrand(TValue (&integrand)(TArgScalar & t)){
		return Integrand<TValue, TArgScalar>(integrand);
	}

	enum Algorithm {
		TRAPEZOIDAL,
		SIMPSON,
		TANH1,
		TANH3,
		TANH5,
		TANH_SINH,
		DEFAULT = SIMPSON
	};

	namespace internal {
		template <typename TValue, int IExponent, typename TArgScalar, typename TFunctor>
		inline TValue tanhFunctor(const TFunctor & f, TArgScalar pos, TArgScalar mean, TArgScalar diffHalf){
			const double posPowExp = IExponent == 1 ? pos : pow(pos, IExponent);
			double ch = cosh(posPowExp);
			if(IExponent == 1)
				return f((TArgScalar)(diffHalf * tanh(posPowExp) + mean)) / ch * ch;
			else
				return f((TArgScalar)(diffHalf * tanh(posPowExp) + mean)) * IExponent * posPowExp / (pos * ch * ch);
		}
	}

	//TODO allow integration according to arbitrary TimePolicy
	template <enum Algorithm EAlgorithm, typename TValue, typename TArgScalar, typename TFunctor>
	inline TValue integrateFunctor(TArgScalar a, TArgScalar b, const TFunctor & f, int numberOfPoints, TValue zero = TValue(0)){
		switch(EAlgorithm){
		case SIMPSON:
			numberOfPoints += (1 - numberOfPoints % 2); // make odd
			break;
		default:
			break;
		}

		TArgScalar diff = b-a, stepSize = diff / (numberOfPoints - 1);
		if(diff == 0) return zero;

		switch(EAlgorithm){
		case TRAPEZOIDAL:
			{
				TValue sum = 0.5 * (f(a) + f(b));

				for(int i = 1; i < numberOfPoints - 1; i++){
					sum += f((TArgScalar)(a + stepSize * i));
				}
				return sum * stepSize;
			}
		case SIMPSON:
			{
				TValue sum = 0.5 * (f(a) + f(b));

				for(int i = 1; i < numberOfPoints - 1; i++){
					sum += f((TArgScalar)(a + stepSize * i)) * (i % 2 == 1 ? 2 : 1);
				}
				return sum * (stepSize * 2. / 3.);
			}
		case TANH1:
		case TANH3:
		case TANH5:
			throw std::runtime_error("TANHX integration scheme not supported yet");
			{
				const double diffHalf = diff /2.;
				const double mean = 0.5 * (a+b);
				const int exponent = 2 * (EAlgorithm - TANH1) + 1;
				TValue sum = internal::tanhFunctor<TValue, exponent>(f, (TArgScalar)0, mean, diffHalf);

				for(int i = numberOfPoints / 2, end = numberOfPoints / 2; i <= end; i++){
					sum += internal::tanhFunctor<TValue, exponent>(f, (TArgScalar) i, mean, diffHalf);
				}
				return sum / diffHalf;
			}
		case TANH_SINH:
			throw std::runtime_error("TANH_SINH integration scheme not supported yet");
			break;
		}
	}
	template <typename TValue, typename TArgScalar, typename TFunctor>
	inline TValue integrateFunctor(TArgScalar a, TArgScalar b, const TFunctor & f, int numberOfPoints, TValue zero = TValue(0)){
		return integrateFunctor<DEFAULT> (a, b, f, numberOfPoints, zero);
	}

	template <enum Algorithm EAlgorithm, typename TValue, typename TArgScalar>
	inline TValue integrateFunction(TArgScalar a, TArgScalar b, TValue (& integrand)(const TArgScalar & t), int numberOfPoints, TValue zero = TValue(0)){
		return integrateFunctor(a, b, createIntegrand(integrand), numberOfPoints, zero);
	}
	template <typename TValue, typename TArgScalar>
	inline TValue integrateFunction(TArgScalar a, TArgScalar b, TValue (& integrand)(const TArgScalar & t), int numberOfPoints, TValue zero = TValue(0)){
		return integrateFunction<DEFAULT>(a, b, integrand, numberOfPoints, zero);
	}
}

#endif /* NUMERICINTEGRATOR_HPP_ */
