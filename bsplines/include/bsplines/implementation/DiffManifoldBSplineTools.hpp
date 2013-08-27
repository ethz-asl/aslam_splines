/*
 * DiffManifoldBSplineTools.hpp
 *
 *  Created on: Aug 27, 2013
 *      Author: hannes
 */

#ifndef DIFFMANIFOLDBSPLINETOOLS_HPP_
#define DIFFMANIFOLDBSPLINETOOLS_HPP_


namespace bsplines {
namespace internal{

template<typename Iterator>
inline void moveIterator(Iterator & it, const Iterator & limit, int steps)
{
	if(steps == 0)
		return;
	if(steps > 0){
		for (int c = steps; c>0; c--){
			if(it == limit)
				break;
			it++;
		}
	}
	else{
		for (int c = -steps; c>0; c--){
			if(it == limit)
				break;
			it--;
		}
	}
}

template<typename Iterator>
inline Iterator getMovedIterator(const Iterator & it, const Iterator & limit, int steps)
{
	Iterator nIt(it);
	moveIterator(nIt, limit, steps);
	return nIt;
}

}
}

#endif /* DIFFMANIFOLDBSPLINETOOLS_HPP_ */
