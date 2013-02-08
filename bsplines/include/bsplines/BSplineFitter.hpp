/*
 * DiffManifoldBSplineFitter.hpp
 *
 *  Created on: May 10, 2012
 *      Author: hannes
 */

#ifndef RIEMANNIANBSPLINEFITTER_HPP_
#define RIEMANNIANBSPLINEFITTER_HPP_

#include "DynamicOrTemplateInt.hpp"
#include "DiffManifoldBSpline.hpp"
#include <sparse_block_matrix/sparse_block_matrix.h>

namespace bsplines{
template <typename TSpline>
class BSplineFitter {
	typedef typename TSpline::point_t point_t;
public :
	SM_DEFINE_EXCEPTION(Exception, std::runtime_error);

	static void initUniformSpline(TSpline & spline, const std::vector<typename TSpline::time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda);
	static void initUniformSplineDense(TSpline & spline, const std::vector<typename TSpline::time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda);
	static void initUniformSplineSparse(TSpline & spline, const std::vector<typename TSpline::time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda);

private:
	static void addCurveQuadraticIntegralDiagTo(const TSpline & spline, const point_t & Wdiag, int derivativeOrder, Eigen::MatrixXd & toMatrix);
	static void addCurveQuadraticIntegralDiagToSparse(const TSpline & spline, const point_t & Wdiag, int derivativeOrder, sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> & toMatrix);

	template<typename M_T>
	static void addOrSetSegmentQuadraticIntegralDiag(const TSpline & spline, const point_t & Wdiag, typename TSpline::SegmentConstIterator segmentIt, int derivativeOrder, M_T toMatrix, bool add);
	template <typename M_T>

	static void computeBijInto(const TSpline & spline, const typename TSpline::SegmentConstIterator & segmentIndex, int columnIndex, M_T B);

	template <typename M_T>
	static void computeMiInto(const TSpline & spline, const typename TSpline::SegmentConstIterator & segmentIndex, M_T & M);
};
}

#include "bsplines/implementation/BSplineFitterImpl.hpp"

#endif /* RIEMANNIANBSPLINEFITTER_HPP_ */
