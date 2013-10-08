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
	template <typename SimpleType_>
	struct SimpleTypeTimePolicy {
		typedef SimpleType_ time_t;
		typedef SimpleType_ duration_t;

		inline static duration_t computeDuration(time_t from, time_t till){
			return till - from;
		}

		inline static duration_t addScaledDuration(time_t from, duration_t dist, int scale){
			return from + dist * scale;
		}

		inline static double divideDurations(duration_t a, duration_t b){
			return (double) a/b;
		}

		inline static time_t linearlyInterpolate(time_t from, time_t till, int segments, int pos)
		{
			if(pos == segments) return till;
			return from + computeDuration(from, till) / segments * pos;
		}

		inline static int getSegmentNumber(time_t from, time_t till, int segments, time_t t)
		{
			double dt = (till - from) / segments;
			return std::floor((t - from) / dt);
		}

		constexpr inline static duration_t getZero(){
			return duration_t(0.0);
		}

		constexpr inline static duration_t getOne(){
			static_assert(!std::is_integral<duration_t>::value, "Please shadow this (getOne) method for your simple type's time policy - preferably with an inline or even better constexpr function.");
			return duration_t(1.0);
		}
	};

} // namespace bsplines

#endif /* SIMPLETYPETIMEPOLICY_HPP_ */
