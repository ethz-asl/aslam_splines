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
#include "KnotArithmetics.hpp"

namespace bsplines{

template <typename TSpline>
class BSplineFitter {
	typedef typename TSpline::point_t point_t;
	typedef typename TSpline::time_t time_t;
	typedef typename TSpline::TimePolicy TimePolicy;
public :
	SM_DEFINE_EXCEPTION(Exception, std::runtime_error);

	enum class Backend {
		DENSE,
		SPARSE,
		DEFAULT = SPARSE
	};

	static void initSpline(TSpline & spline, KnotGenerator<time_t> & knotGenerator, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda, Backend backend = Backend::DEFAULT);
	static void initUniformSpline(TSpline & spline, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda, Backend backend = Backend::DEFAULT);
	static void initUniformSplineDense(TSpline & spline, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda);
	static void initUniformSplineSparse(TSpline & spline, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda);

private:
	static void addCurveQuadraticIntegralDiagTo(const TSpline & spline, const point_t & Wdiag, int derivativeOrder, Eigen::MatrixXd & toMatrix);
	static void addCurveQuadraticIntegralDiagToSparse(const TSpline & spline, const point_t & Wdiag, int derivativeOrder, sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> & toMatrix);

	template<typename M_T>
	static void addOrSetSegmentQuadraticIntegralDiag(const TSpline & spline, const point_t & Wdiag, typename TSpline::SegmentConstIterator segmentIt, int derivativeOrder, M_T toMatrix, bool add);

	template <typename M_T>
	static void computeBijInto(const TSpline & spline, const typename TSpline::SegmentConstIterator & segmentIndex, int columnIndex, M_T B);

	template <typename M_T>
	static void computeMiInto(const TSpline & spline, const typename TSpline::SegmentConstIterator & segmentIndex, M_T & M);

	static void calcFittedControlVerticesDense(TSpline & spline, const KnotIndexResolver<time_t> & knotResolver, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda);
	static void calcFittedControlVerticesSparse(TSpline & spline, const KnotIndexResolver<time_t> & knotResolver, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda);
};

}


#include "bsplines/implementation/BSplineFitterImpl.hpp"

#endif /* RIEMANNIANBSPLINEFITTER_HPP_ */
