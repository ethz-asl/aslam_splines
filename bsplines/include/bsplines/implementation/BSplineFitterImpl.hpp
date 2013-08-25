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

	//TODO improve : support different containers of times and points
	_TEMPLATE
	inline void _CLASS::initUniformSpline(TSpline & spline, const std::vector<time_t> & times, const std::vector<point_t> & points, int numSegments, double lambda, FittingBackend backend){
		const size_t numInterpolationPoints = points.size();
		IntervalUniformKnotGenerator<TimePolicy> knotGenerator(spline.getSplineOrder(), times[0], times[numInterpolationPoints - 1], numSegments);
		initSpline(spline, knotGenerator, times, points, numSegments, lambda, backend);
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

	template<typename TTime>
	struct MapKnotIndexResolver : public KnotIndexResolver<TTime>{
		std::map<TTime, int> knots;
		inline int getKnotIndexAtTime(TTime t) const { return knots.lower_bound(t)->second; };
	};
}

	_TEMPLATE
	void _CLASS::initSpline(TSpline & spline, KnotGenerator<time_t> & knotGenerator, const std::vector<time_t> & times, const std::vector<point_t> & points, int numSegments, double lambda, FittingBackend fittingBackend){
		const int splineOrder = spline.getSplineOrder();
		const size_t numInterpolationPoints = points.size();

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

		if(knotGenerator.hasKnotResolver()){
			for(size_t i = 0; i < K; i ++){
				spline.addKnot(knotGenerator.getNextKnot());
			}
			resolverPtr.init(&knotGenerator.getKnotResolver(), false);
		}else{
			auto mapResolver = new internal::MapKnotIndexResolver<time_t>();
			resolverPtr.init(mapResolver, true);
			for(size_t i = 0; i < K; i ++){
				time_t nextKnot = knotGenerator.getNextKnot();
				spline.addKnot(nextKnot);
				mapResolver->knots[nextKnot] = i;
			}
		}

		spline.init();

		calcFittedControlVertices(spline, *resolverPtr, times, points, lambda, fittingBackend);
	}

	_TEMPLATE
	void _CLASS::fitSpline(TSpline & spline, const std::vector<time_t> & times, const std::vector<point_t> & points, double lambda, FittingBackend fittingBackend){
		spline.assertEvaluable();

		internal::MapKnotIndexResolver<time_t> mapResolver;
		int i = 0;
		for(auto it = spline.getFirstRelevantSegmentByLast(spline.getSegmentIterator(times[0])), end = ++spline.getSegmentIterator(times[times.size() - 1]); it != end; ++it){
			time_t nextKnot = it->getKnot();
			spline.addKnot(nextKnot);
			mapResolver.knots[nextKnot] = i++;
		}
		calcFittedControlVertices(spline, mapResolver, times, points, lambda, fittingBackend);
	}

	_TEMPLATE
	void _CLASS::initUniformSplineDense(TSpline & spline, const std::vector<time_t> & times, const std::vector<point_t> & points, int numSegments, double lambda)
	{
		initUniformSpline(spline, times, points, numSegments, lambda, FittingBackend::DENSE);
	}

	_TEMPLATE
	void _CLASS::initUniformSplineSparse(TSpline & spline, const std::vector<time_t> & times, const std::vector<point_t> & points, int numSegments, double lambda)
	{
		initUniformSpline(spline, times, points, numSegments, lambda, FittingBackend::SPARSE);
	}

	_TEMPLATE
	void _CLASS::calcFittedControlVertices(TSpline & spline, const KnotIndexResolver<time_t> & knotResolver, const std::vector<time_t> & times, const std::vector<point_t> & points, double lambda, FittingBackend fittingBackend, int fixNFirstRelevantControlVertices){
		switch(fittingBackend){
			case FittingBackend::DENSE:
				calcFittedControlVertices<FittingBackend::DENSE>(spline, knotResolver, times, points, lambda);
				break;
			case FittingBackend::SPARSE:
				calcFittedControlVertices<FittingBackend::SPARSE>(spline, knotResolver, times, points, lambda);
				break;
		}
	}

	namespace internal {
		template< enum FittingBackend FittingBackend_> struct FittingBackendTraits {
		};

		template<> struct FittingBackendTraits<FittingBackend::DENSE> {
			typedef Eigen::MatrixXd TypeA;
			typedef Eigen::VectorXd TypeB;

			inline TypeA createA(unsigned int constraintSize, unsigned int coefficientDim, unsigned int D) {
				return TypeA::Zero(constraintSize, coefficientDim);
			}
			inline TypeB createB(int constraintSize) {
				return TypeB::Zero(constraintSize);
			}

			inline auto blockA(TypeA & A, int row, int col, int D) -> decltype(A.block(0, 0, 0, 0)){
				return A.block(row * D, col * D, D, D);
			}
			inline typename TypeB::SegmentReturnType segmentB(TypeB & b, int row, int D){
				return b.segment(row * D, D);
			}

			inline Eigen::VectorXd solve(TypeA & A, TypeB & b){
				return A.ldlt().solve(b);
			}
		};

		template<> struct FittingBackendTraits<FittingBackend::SPARSE> {
			typedef sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> TypeA;
			typedef sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> TypeB;

			std::vector<int> rows;
			std::vector<int> cols;

			const bool allocateBlock = true;

			inline TypeA createA(unsigned int constraintSize, unsigned int coefficientDim, unsigned int D) {
				for(unsigned int i = D; i <= constraintSize; i+=D) rows.push_back(i);
				for(unsigned int i = D; i <= coefficientDim; i+=D) cols.push_back(i);
				return TypeA(rows,cols, true);
			}

			inline TypeB createB(int constraintSize) {
				std::vector<int> bcols(1);
				bcols[0] = 1;
				return TypeB(rows,bcols, true);
			}

			inline typename TypeA::SparseMatrixBlock & blockA(TypeA & A, int row, int col, int D){
				return *A.block(row, col, allocateBlock );
			}
			inline typename TypeB::SparseMatrixBlock & segmentB(TypeB & b, int row, int D){
				return *b.block(row, 0, allocateBlock);
			}

			inline Eigen::VectorXd solve(TypeA & A, TypeB & b){
				// solve:
				sparse_block_matrix::LinearSolverCholmod<Eigen::MatrixXd> solver;
				solver.init();

				Eigen::VectorXd c(A.rows());
				c.setZero();
				Eigen::VectorXd b_dense = b.toDense();
				bool result = solver.solve(A,&c[0],&b_dense[0]);
				if(!result) {
					c.setZero();
					// fallback => use nonsparse solver:
					std::cout << "Fallback to Dense Solver" << std::endl;
					Eigen::MatrixXd Adense = A.toDense();
					c = Adense.ldlt().solve(b_dense);
				}
				return c;
			}
		};
	}

	_TEMPLATE
	template <enum FittingBackend FittingBackend_>
	void _CLASS::calcFittedControlVertices(TSpline & spline, const KnotIndexResolver<time_t> & knotResolver, const std::vector<time_t> & times, const std::vector<point_t> & points, double lambda, int fixNFirstRelevantControlVertices)
	{
		const int splineOrder = spline.getSplineOrder();
		const int minKnotIndex = knotResolver.getKnotIndexAtTime(times[0]) - (splineOrder - 1); // get minimal relevant knot's index
		const int controlVertexIndexOffset = minKnotIndex + fixNFirstRelevantControlVertices;

		const typename TSpline::point_t* fixControlVertices[fixNFirstRelevantControlVertices];
		if(fixNFirstRelevantControlVertices){
			auto it = spline.getFirstRelevantSegmentByLast(spline.getSegmentIterator(times[0]));
			for(int i = 0; i < fixNFirstRelevantControlVertices; ++i){
				fixControlVertices[i] = &it->getControlVertex();
				++it;
			}
		}

		const size_t numInterpolationPoints = points.size();

		// What is the vector coefficient dimension
		const size_t D = spline.getPointSize();

		// How many coefficients are required for one time segment?
		const size_t C = spline.getNumControlVertices() - controlVertexIndexOffset; // we are fitting the tail of the spline

		// Now we have to solve an Ax = b linear system to determine the correct coefficient vectors.
		size_t coefficientDim = C * D;

		const int constraintSize = numInterpolationPoints * D;

		internal::FittingBackendTraits<FittingBackend_> backend;

		auto A = backend.createA(constraintSize, coefficientDim, D);
		auto b = backend.createB(constraintSize);

		int brow = 0;
		// Add the position constraints.
		for(size_t i = 0; i < numInterpolationPoints; i++)
		{
			int knotIndex = knotResolver.getKnotIndexAtTime(times[i]) - controlVertexIndexOffset;
			//TODO optimize : a uniform time spline evaluator could be much faster here!
			typename TSpline::SplineOrderVector bi = spline.template getEvaluatorAt<1>(times[i]).getLocalBi();

			knotIndex -=  splineOrder - 1; // get first relevant knot's index
			if(i == numInterpolationPoints - 1){ //TODO improve: remove this irregularity
				knotIndex --;
			}

			for(int j = 0; j < splineOrder; j ++){
				int col = knotIndex + j;
				if(col >= 0){ // this control vertex is not fixed
					if(bi[j] != 0.0)
						backend.blockA(A, brow, col, D).diagonal().setConstant(bi[j]);
				}
				else{
					backend.segmentB(b, brow, D) += (*fixControlVertices[fixNFirstRelevantControlVertices + j]) * bi[j];
				}
			}

			if(knotIndex >= 0)
				backend.segmentB(b, brow, D) = points[i];
			else{
				backend.segmentB(b, brow, D) += points[i];
			}
			++brow;
		}

		auto At = A.transpose();
		b = (At * b).eval();
		A = (At * A).eval();

		if(lambda != 0.0){
			// Add the motion constraint.
			point_t W = point_t::Constant(D, lambda);
			addCurveQuadraticIntegralDiagTo(spline, W, 2, A);
		}

		// Solve for the coefficient vector.
		Eigen::VectorXd c = backend.solve(A, b);

		spline.setControlVertices(c, times[0]);
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
	void _CLASS::addCurveQuadraticIntegralDiagTo(const TSpline & spline, const point_t & Wdiag, int derivativeOrder, sparse_block_matrix::SparseBlockMatrix<Eigen::MatrixXd> & toMatrix)
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
