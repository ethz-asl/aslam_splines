/*
 * SimpleTypeTimePolicy.hpp
 *
 *  Created on: May 10, 2012
 *      Author: hannes
 */

#ifndef SIMPLETYPETIMEPOLICY_HPP_
#define SIMPLETYPETIMEPOLICY_HPP_

#include <cmath>

namespace bsplines {
	template <typename SIMPLE_TYPE>
	struct SimpleTypeTimePolicy {
		typedef SIMPLE_TYPE time_t;
		typedef SIMPLE_TYPE duration_t;

		inline static duration_t computeDuration(time_t from, time_t till){
			return till - from;
		}

		inline static double divideDurations(duration_t a, duration_t b){
			return (double) a/b;
		}

		inline static time_t linearlyInterpolate(time_t from, time_t till, int segments, int pos)
		{
			return from + computeDuration(from, till) / segments * pos;
		}

		inline static int getSegmentNumber(time_t from, time_t till, int segments, time_t t)
		{
			double dt = (till - from) / segments;
			return std::floor((t - from) / dt);
		}

		inline static duration_t getZero(){
			return 0.0;
		}

		inline static duration_t getOne(){
			return 1.0;
		}
	};
}

#endif /* SIMPLETYPETIMEPOLICY_HPP_ */
