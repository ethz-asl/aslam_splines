/*
 * DiffManifoldBSplineFitter.hpp
 *
 *  Created on: May 10, 2012
 *      Author: hannes
 */

#ifndef RIEMANNIANBSPLINEFITTERIMPL_HPP_
#define RIEMANNIANBSPLINEFITTERIMPL_HPP_

#include "../DynamicOrTemplateInt.hpp"
#include "../DiffManifoldBSpline.hpp"

#include <memory>
#include <map>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/QR>
#include <sparse_block_matrix/linear_solver_cholmod.h>

namespace bsplines{
#define _TEMPLATE template <typename TSpline>
#define _CLASS BSplineFitter <TSpline>

	//TODO improve : support different containers of times and interpolationPoints
	_TEMPLATE
	inline void _CLASS::initUniformSpline(TSpline & spline, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda, Backend backend){
		const size_t numInterpolationPoints = interpolationPoints.size();
		IntervalUniformKnotGenerator<TimePolicy> knotGenerator(spline.getSplineOrder(), times[0], times[numInterpolationPoints - 1], numSegments);
		initSpline(spline, knotGenerator, times, interpolationPoints, numSegments, lambda, backend);
	}

namespace internal{
	template<typename T>
	class PointerGuard {
		bool _deleteIt;
		T * _ptr;
	 public:
		PointerGuard() : _deleteIt(false), _ptr(nullptr) {}
		inline void init(T* ptr, bool deleteIt = true) {
			_deleteIt = deleteIt;
			_ptr = ptr;
		}
		inline ~PointerGuard(){
			if(_deleteIt) delete _ptr;
		}
		inline T& operator *(){
			return *_ptr;
		}
	};
}

	_TEMPLATE
	void _CLASS::initSpline(TSpline & spline, KnotGenerator<time_t> & knotGenerator, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda, Backend backend){
		const int splineOrder = spline.getSplineOrder();
		const size_t numInterpolationPoints = interpolationPoints.size();

		SM_ASSERT_EQ(Exception, times.size(), numInterpolationPoints, "The number of times and the number of interpolation points must be equal");
		SM_ASSERT_GE(Exception, numInterpolationPoints, 2, "There must be at least two times/points");
		SM_ASSERT_GE(Exception, numSegments, 1, "There must be at least one time segment");
		for(size_t i = 1; i < numInterpolationPoints; i++)
		{
			SM_ASSERT_LE(Exception, times[i-1], times[i],
					"The time sequence must be nondecreasing. time " << i
					<< " was not greater than or equal to time " << (i-1));
		}

		// How many knots are required for one time segment?
		const size_t K = knot_arithmetics::getNumKnotsRequired(numSegments, splineOrder);

		internal::PointerGuard<const KnotIndexResolver<time_t> > resolverPtr;

		struct MapKnotResolver : public KnotIndexResolver<time_t>{
			std::map<time_t, int> knots;
			inline int getKnotIndexAtTime(time_t t) const { return knots.lower_bound(t)->second; };
		};

		if(knotGenerator.hasKnotResolver()){
			for(size_t i = 0; i < K; i ++){
				spline.addKnot(knotGenerator.getNextKnot());
			}
			resolverPtr.init(&knotGenerator.getKnotResolver(), false);
		}else{
			MapKnotResolver * mapResolver = new MapKnotResolver();
			resolverPtr.init(mapResolver, true);
			for(size_t i = 0; i < K; i ++){
				time_t nextKnot = knotGenerator.getNextKnot();
				spline.addKnot(nextKnot);
				mapResolver->knots[nextKnot] = i;
			}
		}

		spline.init();

		switch(backend){
			case Backend::SPARSE:
				calcFittedControlVerticesSparse(spline, *resolverPtr, times, interpolationPoints, numSegments, lambda);
				break;
			case Backend::DENSE:
				calcFittedControlVerticesDense(spline, *resolverPtr, times, interpolationPoints, numSegments, lambda);
				break;
		}
	}

	_TEMPLATE
	void _CLASS::initUniformSplineDense(TSpline & spline, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda)
	{
		initUniformSpline(spline, times, interpolationPoints, numSegments, lambda, Backend::DENSE);
	}

	_TEMPLATE
	void _CLASS::initUniformSplineSparse(TSpline & spline, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda)
	{
		initUniformSpline(spline, times, interpolationPoints, numSegments, lambda, Backend::SPARSE);
	}

	_TEMPLATE
	void _CLASS::calcFittedControlVerticesDense(TSpline & spline, const KnotIndexResolver<time_t> & knotResolver, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda)
	{
		const int splineOrder = spline.getSplineOrder();
		const size_t numInterpolationPoints = interpolationPoints.size();

		// What is the vector coefficient dimension
		const size_t D = spline.getPointSize();
		// How many coefficients are required for one time segment?
		const size_t C = knot_arithmetics::getNumControlVerticesRequired(numSegments, splineOrder);

		// Now we have to solve an Ax = b linear system to determine the correct coefficient vectors.
		size_t coefficientDim = C * D;

		const int numConstraints = numInterpolationPoints;
		const int constraintSize = numConstraints * D;

		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(constraintSize, coefficientDim);
		Eigen::VectorXd b = Eigen::VectorXd::Zero(constraintSize);

		int brow = 0;
		// Add the position constraints.
		for(size_t i = 0; i < numInterpolationPoints; i++)
		{
			int knotIndex = knotResolver.getKnotIndexAtTime(times[i]);
			//TODO optimize : a uniform time spline evaluator could be much faster here!
			typename TSpline::SplineOrderVector bi = spline.template getEvaluatorAt<1>(times[i]).getLocalBi();

			knotIndex -=  splineOrder - 1;
			if(i == numInterpolationPoints - 1){ //TODO improve: remove this irregularity
				knotIndex --;
			}

			knotIndex *= D;
			for(int j = 0; j < splineOrder; j ++){
				if(bi[j] != 0.0)
					A.block(brow, knotIndex + j * D, D, D).diagonal().setConstant(bi[j]);
			}

			b.segment(brow,D) = interpolationPoints[i];
			brow += D;
		}

		b = (A.transpose() * b).eval();
		A = (A.transpose() * A).eval();

		if(lambda != 0.0){
			// Add the motion constraint.
			point_t W = point_t::Constant(D, lambda);
			addCurveQuadraticIntegralDiagTo(spline, W, 2, A);
		}

		// Solve for the coefficient vector.
		Eigen::VectorXd c = A.ldlt().solve(b);

		size_t i = 0;
		for(typename TSpline::SegmentIterator it = spline.getAbsoluteBegin(), end = spline.getAbsoluteEnd(); it != end && i < coefficientDim; it ++){
			it->getControlVertex() = c.segment(i, D);
			spline.getManifold().projectIntoManifold(it->getControlVertex());
			i+= D;
		}
	}

	_TEMPLATE
	void _CLASS::calcFittedControlVerticesSparse(TSpline & spline, const KnotIndexResolver<time_t> & knotResolver, const std::vector<time_t> & times, const std::vector<point_t> & interpolationPoints, int numSegments, double lambda)
	{
		const int splineOrder = spline.getSplineOrder();
		const size_t numInterpolationPoints = interpolationPoints.size();

		// What is the vector coefficient dimension
		const size_t D = spline.getPointSize();
		// How many coefficients are required for one time segment?
		const size_t C = knot_arithmetics::getNumControlVerticesRequired(numSegments, splineOrder);


		// Initialize a uniform knot sequence
		const knot_arithmetics::UniformTimeCalculator<typename TSpline::TimePolicy> timeCalculator(splineOrder, times[0], times[numInterpolationPoints - 1], numSegments);

		// Now we have to solve an Ax = b linear system to determine the correct coefficient vectors.
		size_t coefficientDim = C * D;

		// define the structure:
		std::vector<int> rows;
		std::vector<int> cols;

		for (unsigned int i = 1; i <= numInterpolationPoints; i++) rows.push_back(i*D);
		for(unsigned int i = 1; i <= C; i++) cols.push_back(i*D);

		std::vector<int> bcols(1);
		bcols[0] = 1;

		sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> A(rows,cols, true);
		sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> b(rows,bcols, true);

		int brow = 0;
		// try to fill the matrix:
		for(unsigned int i = 0; i < numInterpolationPoints; i++) {
			int knotIndex = knotResolver.getKnotIndexAtTime(times[i]);
			//TODO optimize : a uniform time spline evaluator could be much faster here!
			typename TSpline::SplineOrderVector bi = spline.template getEvaluatorAt<1>(times[i]).getLocalBi();

			knotIndex -=  splineOrder - 1;
			if(i == numInterpolationPoints - 1){ //TODO improve: remove this irregularity
				knotIndex --;
			}

			const bool allocateBlock = true;

			knotIndex *= D;

			// the n'th order spline needs n column blocks (n*D columns)
			for(int j = 0; j < splineOrder; j++) {
				Eigen::MatrixXd & Ai = *A.block(brow/D,knotIndex/D + j,allocateBlock );
				Ai.diagonal().setConstant(bi[j]);
			}

			*b.block(brow/D,0,allocateBlock ) = interpolationPoints[i];

			brow += D;
		}

		sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> At(cols,rows, true);
		sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> * Atp = &At;
		A.transpose(Atp);

		// A'b
		sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> Ab(cols,bcols, true);
		sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> * Abp = &Ab;
		Atp->multiply(Abp, &b);

		// A'A
		sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> AtA(cols,cols, true);
		sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> * AtAp = &AtA;
		Atp->multiply(AtAp, &A);

		if(lambda != 0.0){
			// Add the motion constraint.
			point_t W = point_t::Constant(D, lambda);
			addCurveQuadraticIntegralDiagToSparse(spline, W, 2, *AtAp);
		}

		// solve:
		sparse_block_matrix::LinearSolverCholmod<Eigen::MatrixXd> solver;
		solver.init();

		Eigen::VectorXd c(AtAp->rows());
		c.setZero();
		Eigen::VectorXd b_dense = Abp->toDense();

		bool result = solver.solve(*AtAp,&c[0],&b_dense[0]);
		if(!result) {
			c.setZero();
			// fallback => use nonsparse solver:
					std::cout << "Fallback to Dense Solver" << std::endl;
					Eigen::MatrixXd Adense = AtAp->toDense();
					c = Adense.ldlt().solve(b_dense);
		}

		size_t i = 0;
		for(typename TSpline::SegmentIterator it = spline.getAbsoluteBegin(), end = spline.getAbsoluteEnd(); it != end && i < coefficientDim; it ++){
			it->getControlVertex() = c.segment(i, D);
			spline.getManifold().projectIntoManifold(it->getControlVertex());
			i+= D;
		}
	}



	_TEMPLATE
	void _CLASS::addCurveQuadraticIntegralDiagTo(const TSpline & spline, const point_t & Wdiag, int derivativeOrder, Eigen::MatrixXd & toMatrix)
	{
		SM_ASSERT_EQ(Exception,Wdiag.rows(), spline.getPointSize(), "Wdiag must be of control point size");

		int qiSize = spline.getSplineOrder() * spline.getPointSize();
		int brow = 0;
		for(typename TSpline::SegmentConstIterator sIt = spline.begin(), end = spline.end(); sIt != end; sIt++)
		{
			addOrSetSegmentQuadraticIntegralDiag(spline, Wdiag, sIt, derivativeOrder, toMatrix.block(brow, brow, qiSize, qiSize), true);
			brow += spline.getPointSize();
		}
	}

	// sparse curveQuaddraticIntegral:
	_TEMPLATE
	void _CLASS::addCurveQuadraticIntegralDiagToSparse(const TSpline & spline, const point_t & Wdiag, int derivativeOrder, sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> & toMatrix)
	{
		SM_ASSERT_EQ(Exception,Wdiag.rows(), spline.getPointSize(), "Wdiag must be of control point size");

		BOOST_AUTO(qiSize, spline.getSplineOrder() * spline.getPointSize());
		typedef BOOST_TYPEOF_TPL(qiSize) QiSize;
		typedef Eigen::Matrix<double, QiSize::VALUE, QiSize::VALUE> Q_T;
		int brow = 0;
		Q_T Q((int)qiSize, (int)qiSize);

		int D = spline.getPointSize();
		int blocksInBlock = spline.getSplineOrder();
//
		const bool allocateBlock = true;
		int s = 0;
		for(typename TSpline::SegmentConstIterator sIt = spline.begin(), end = spline.end(); sIt != end; sIt++)
		{
			addOrSetSegmentQuadraticIntegralDiag<Q_T &>(spline, Wdiag, sIt, derivativeOrder, Q, false);
			brow += spline.getPointSize();
			// place the DxD blocks in the blocksInBlock x blocksInBlock blocks:
			for(int i = 0; i < blocksInBlock; i++) {
				for(int j = 0; j < blocksInBlock; j++) {
					Eigen::MatrixXd & Qi = *toMatrix.block(s+i, s+j, allocateBlock);
					Qi += Q.block(i*D,j*D,D,D);
				}
			}
			s++;
		}
	}

	_TEMPLATE
	template<typename M_T>
	inline void _CLASS::addOrSetSegmentQuadraticIntegralDiag(const TSpline & spline, const point_t & Wdiag, typename TSpline::SegmentConstIterator segmentIt, int derivativeOrder, M_T toMatrix, bool add)
	{
		const int D = spline.getPointSize();
		const int splineOrder = spline.getSplineOrder();
		SM_ASSERT_GE_LT(Exception, segmentIt.getKnot(), spline.getMinTime(), spline.getMaxTime(), "Out of range");
		SM_ASSERT_EQ(Exception, Wdiag.rows(), D, "Wdiag must be of control point size");

		typename TSpline::SplineOrderSquareMatrix Dm(splineOrder, splineOrder);
		Dm.setZero(splineOrder, splineOrder);
		spline.computeDiiInto(segmentIt, Dm);
		typename TSpline::SplineOrderSquareMatrix V(splineOrder, splineOrder);
		V.setZero(splineOrder, splineOrder);
		spline.computeViInto(segmentIt, V);

		// Calculate the appropriate derivative version of V
		// using the matrix multiplication version of the derivative.
		for(int i = 0; i < derivativeOrder; i++)
		{
			V = (Dm.transpose() * V * Dm).eval();
		}

		BOOST_AUTO(splineOrderTimesPointSize, spline.getSplineOrder() * spline.getPointSize());
		typedef BOOST_TYPEOF_TPL(splineOrderTimesPointSize) SplineOrderTimesPointSize;
		typedef Eigen::Matrix<double, SplineOrderTimesPointSize::VALUE, SplineOrderTimesPointSize::VALUE> MType;

		MType WV((int)splineOrderTimesPointSize, (int)splineOrderTimesPointSize);

		WV.setZero(splineOrderTimesPointSize, splineOrderTimesPointSize);

		for(int d = 0; d < D; d++)
		{
			WV.block(splineOrder*d, splineOrder*d, splineOrder, splineOrder) = Wdiag(d) * V;
		}

		MType M((int)splineOrderTimesPointSize, (int)splineOrderTimesPointSize);
		computeMiInto(spline, segmentIt, M);

		if(add) toMatrix += M.transpose() * WV * M;
		else toMatrix = M.transpose() * WV * M;
	}

	_TEMPLATE
	template <typename M_T>
	void _CLASS::computeBijInto(const TSpline & spline, const typename TSpline::SegmentConstIterator & segmentIndex, int columnIndex, M_T B)
	{
		const int D = spline.getPointSize();
		const int splineOrder = spline.getSplineOrder();

		for(int i = 0; i < D; i++)
		{
			B.block(i*splineOrder,i,splineOrder,1) = segmentIndex->getBasisMatrix().col(columnIndex);
		}
	}

	_TEMPLATE
	template <typename M_T>
	void _CLASS::computeMiInto(const TSpline & spline, const typename TSpline::SegmentConstIterator & segmentIndex, M_T & M)
	{
		const int D = spline.getPointSize();
		const int splineOrder = spline.getSplineOrder();
		const int splineOrderTimesPointSize = splineOrder * D;
		M.setZero();

		for(int j = 0; j < splineOrder; j++)
		{
			computeBijInto(spline, segmentIndex, j, M.block(0, j*D, splineOrderTimesPointSize, D));
		}
	}
}

#undef _TEMPLATE
#undef _CLASS


#endif /* RIEMANNIANBSPLINEFITTERIMPL_HPP_ */
