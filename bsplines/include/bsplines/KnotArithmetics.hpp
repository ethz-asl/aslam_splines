/*
 * BSplineKnotArithmetics.hpp
 *
 *  Created on: May 10, 2012
 *      Author: hannes
 */

#ifndef BSPLINEKNOTARITHMETICS_HPP_
#define BSPLINEKNOTARITHMETICS_HPP_

#include <sm/assert_macros.hpp>

namespace bsplines{
namespace KnotArithmetics{
inline int getNumControlVerticesRequired(int numTimeSegments, int splineOrder)
{
	return numTimeSegments + splineOrder - 1;
}

inline int getNumKnotsRequired(int numTimeSegments, int splineOrder)
{
	return getNumControlVerticesRequired(numTimeSegments, splineOrder) + splineOrder;
}

inline int getNumValidTimeSegments(int numKnots, int splineOrder)
{
	int nv = numKnots - 2 * splineOrder + 1;
	return std::max(nv,0);
}

template <typename TTimePolicy>
class UniformTimeCalculator {
public:
	typedef typename TTimePolicy::time_t time_t;
private:
	const time_t _minValidTime, _maxValidTime;
	const int _numSegments;
	const int _splineOrder;
public:
	inline UniformTimeCalculator(int splineOrder, const time_t minValidTime, const time_t maxValidTime, int numSegments) :
		_minValidTime(minValidTime), _maxValidTime(maxValidTime),
		_numSegments(numSegments), _splineOrder(splineOrder)
	{
	}

	time_t getTimeByKnotIndex(int knotIndex) const {
		SM_ASSERT_GE_LT_DBG(std::runtime_error, knotIndex, 0, getNumKnotsRequired(_numSegments, _splineOrder), " knot indices are only valid in the given range");
		return TTimePolicy::linearlyInterpolate(_minValidTime, _maxValidTime, _numSegments, knotIndex - _splineOrder + 1);
	}

	int getKnotIndexAtTime(time_t t) const {
		return TTimePolicy::getSegmentNumber(_minValidTime, _maxValidTime, _numSegments, t) + _splineOrder - 1;
	}
};
}
}


#endif /* BSPLINEKNOTARITHMETICS_HPP_ */
