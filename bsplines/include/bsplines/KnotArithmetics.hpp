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
	virtual bool hasNextKnot() const = 0;
	virtual Time_ getNextKnot() = 0;
	virtual void jumpOverNextKnots(int amount) { for (int i = 0; i < amount; i++) getNextKnot(); };
	virtual bool hasKnotIndexResolver() { return false; }
	virtual const KnotIndexResolver<Time_>& getKnotIndexResolver() const { throw std::runtime_error("unsupported operation");};
	virtual bool supportsAppending() const { return false; }
	virtual void extendBeyondTime(const Time_ beyondThisTime) { throw std::runtime_error("unsupported operation"); }
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
	inline IntervalUniformKnotGenerator(int splineOrder, time_t tMin, time_t tMax, int numSegments) : _timeCalculator(splineOrder, tMin, tMax, numSegments), _nextIndex(0), _numKnots(knot_arithmetics::getNumKnotsRequired(numSegments, splineOrder)) {
		SM_ASSERT_GT(std::runtime_error, tMax, tMin, "The max time is less than the min time");
	}
	virtual bool hasNextKnot() const { return _nextIndex < _numKnots; };
	virtual time_t getNextKnot() { return _timeCalculator.getTimeByKnotIndex(_nextIndex++); };
	virtual bool hasKnotIndexResolver() { return true; }
	virtual const KnotIndexResolver<time_t>& getKnotIndexResolver() const { return _timeCalculator; };
	virtual ~IntervalUniformKnotGenerator() {}
 private:
	const knot_arithmetics::UniformTimeCalculator<TimePolicy_> _timeCalculator;
	int _nextIndex;
	const int _numKnots;
};

template <typename TimePolicy_>
class DeltaUniformKnotGenerator : public KnotGenerator<typename TimePolicy_::time_t> {
 public:
	typedef typename TimePolicy_::time_t time_t;
	typedef typename TimePolicy_::duration_t duration_t;
	inline DeltaUniformKnotGenerator(int splineOrder, time_t start, time_t beyondThisTime, duration_t delta, bool splineIsAlreadyInitialized = false)
		: _splineOrder(splineOrder), _delta(delta), _lastKnot(TimePolicy_::addScaledDuration(start, delta, splineIsAlreadyInitialized ? -1: -(knot_arithmetics::getNumRequiredPreambleKnots(splineOrder) + 1))), _beyondThisTime(getExtendedEndTime(beyondThisTime)) {}
	DeltaUniformKnotGenerator() = default;
	virtual bool hasNextKnot() const { return _lastKnot < _beyondThisTime; };
	virtual time_t getNextKnot() { return (_lastKnot = TimePolicy_::addScaledDuration(_lastKnot, _delta, 1)); };
	virtual void extendBeyondTime(const time_t beyondThisTime) { _beyondThisTime = getExtendedEndTime(beyondThisTime); }
	virtual bool supportsAppending() const { return true; }
	virtual ~DeltaUniformKnotGenerator() {}
 private:
	inline time_t getExtendedEndTime(time_t time) { return TimePolicy_::addScaledDuration(time, _delta, _splineOrder - 1); }
	int _splineOrder;
	duration_t _delta;
	time_t _lastKnot, _beyondThisTime;
};

}


#endif /* BSPLINEKNOTARITHMETICS_HPP_ */
