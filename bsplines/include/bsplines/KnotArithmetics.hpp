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

template <typename Time_>
class KnotIndexResolver {
 public:
	virtual int getKnotIndexAtTime(Time_ t) const = 0;
	virtual ~KnotIndexResolver() {}
};

template <typename Time_>
class KnotGenerator{
 public:
	virtual Time_ getNextKnot() = 0;
	virtual bool hasKnotResolver() { return false; }
	virtual const KnotIndexResolver<Time_>& getKnotResolver() const { throw std::runtime_error("unsupported operation");};
	virtual bool supportsAppending() const { return false; }
	virtual ~KnotGenerator() {}
};

namespace knot_arithmetics{

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

inline int getNumRequiredPreambleKnots(int splineOrder){
	return splineOrder - 1;
}

template <typename TTimePolicy>
class UniformTimeCalculator : public KnotIndexResolver<typename TTimePolicy::time_t>{
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
		return TTimePolicy::linearlyInterpolate(_minValidTime, _maxValidTime, _numSegments, knotIndex - getNumRequiredPreambleKnots(_splineOrder));
	}

	virtual int getKnotIndexAtTime(time_t t) const {
		return TTimePolicy::getSegmentNumber(_minValidTime, _maxValidTime, _numSegments, t) + getNumRequiredPreambleKnots(_splineOrder);
	}

	virtual ~UniformTimeCalculator(){}
};

}

template <typename TimePolicy_>
class IntervalUniformKnotGenerator : public KnotGenerator<typename TimePolicy_::time_t> {
 public:
	typedef typename TimePolicy_::time_t time_t;
	inline IntervalUniformKnotGenerator(int splineOrder, time_t tMin, time_t tMax, int numSegments) : _timeCalculator(splineOrder, tMin, tMax, numSegments), _nextIndex(0) {
		SM_ASSERT_GT(std::runtime_error, tMax, tMin, "The max time is less than the min time");
	}
	virtual time_t getNextKnot() { return _timeCalculator.getTimeByKnotIndex(_nextIndex++); };
	virtual bool hasKnotResolver() { return true; }
	virtual const KnotIndexResolver<time_t>& getKnotResolver() const { return _timeCalculator; };
	virtual ~IntervalUniformKnotGenerator() {}
 private:
	knot_arithmetics::UniformTimeCalculator<TimePolicy_> _timeCalculator;
	int _nextIndex;
};

template <typename TimePolicy_>
class DeltaUniformKnotGenerator : public KnotGenerator<typename TimePolicy_::time_t> {
 public:
	typedef typename TimePolicy_::time_t time_t;
	typedef typename TimePolicy_::duration_t duration_t;
	inline DeltaUniformKnotGenerator(time_t start, duration_t delta, int splineOrder, bool splineIsAlreadyInitialized = false) : _lastKnot(splineIsAlreadyInitialized ? start : TimePolicy_::addScaledDuration( start, delta, -(knot_arithmetics::getNumRequiredPreambleKnots(splineOrder) + 1))), _delta(delta) {}
	virtual time_t getNextKnot() { return (_lastKnot = TimePolicy_::addScaledDuration(_lastKnot, _delta, 1)); };
	virtual bool supportsAppending() const { return true; }
	virtual ~DeltaUniformKnotGenerator() {}
 private:
	time_t _lastKnot;
	duration_t _delta;
};

}


#endif /* BSPLINEKNOTARITHMETICS_HPP_ */
