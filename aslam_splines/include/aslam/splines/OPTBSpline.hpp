/*
 * OPTBSpline.hpp
 *
 *  Created on: 05.08.2012
 *      Author: hannes
 */

#ifndef OPTBSPLINE_HPP_
#define OPTBSPLINE_HPP_

#include <aslam/backend/DesignVariableMappedVector.hpp>
#include <bsplines/DiffManifoldBSpline.hpp>
#include <Eigen/StdVector>
#include <boost/ptr_container/ptr_vector.hpp>
#include <aslam/backend/TransformationExpression.hpp>
#include <aslam/backend/RotationExpression.hpp>
#include <aslam/backend/EuclideanExpression.hpp>
#include <aslam/backend/VectorExpression.hpp>
#include <aslam/backend/VectorExpressionNode.hpp>

#include <bsplines/manifolds/LieGroup.hpp>


namespace aslam {
	namespace splines {
		namespace internal {
			template <typename TDiffManifoldBSplineConfiguration>
			struct DesignVariableSegmentBSplineConf : public TDiffManifoldBSplineConfiguration {

				typedef TDiffManifoldBSplineConfiguration ParentConf;
				typedef DesignVariableSegmentBSplineConf<TDiffManifoldBSplineConfiguration> Conf;

				typedef typename TDiffManifoldBSplineConfiguration::ManifoldConf ManifoldConf;
				typedef ::bsplines::DiffManifoldBSpline<TDiffManifoldBSplineConfiguration, Conf> BSpline;

				DesignVariableSegmentBSplineConf(TDiffManifoldBSplineConfiguration parentConf) : ParentConf(parentConf) {}
			};
		}
	}
}

namespace bsplines {
	namespace internal {
		template <typename TDiffManifoldBSplineConfiguration>
		struct SegmentData< ::aslam::splines::internal::DesignVariableSegmentBSplineConf<TDiffManifoldBSplineConfiguration> > : public SegmentData<TDiffManifoldBSplineConfiguration>, aslam::backend::DesignVariable{
		private:
			typedef SegmentData<TDiffManifoldBSplineConfiguration> parent_t;
			typedef SegmentData< ::aslam::splines::internal::DesignVariableSegmentBSplineConf<TDiffManifoldBSplineConfiguration> > this_t;

		public:
			typedef typename TDiffManifoldBSplineConfiguration::Manifold::point_t point_t;
			typedef typename TDiffManifoldBSplineConfiguration::Manifold::tangent_vector_t tangent_vector_t;

		private:
			point_t _p_v;

			typename TDiffManifoldBSplineConfiguration::Dimension dimension;

			enum {
				Dimension = TDiffManifoldBSplineConfiguration::Dimension::VALUE
			};

		public:
			SegmentData(const TDiffManifoldBSplineConfiguration & configuration, const time_t & time, const point_t &point) : parent_t(configuration, time, point), dimension(configuration.getDimension()) {}

			virtual int minimalDimensionsImplementation() const { return dimension.getValue(); };

			/// \brief Update the design variable.
			virtual void updateImplementation(const double * dp, int size){
				if(!TDiffManifoldBSplineConfiguration::Dimension::IS_DYNAMIC)
					SM_ASSERT_EQ_DBG(std::runtime_error, Dimension, size, "");
				_p_v = this->getControlVertex();
				Eigen::Map<const tangent_vector_t> dpV(dp, size);
				manifolds::internal::DiffManifoldPointUpdateTraits<typename TDiffManifoldBSplineConfiguration::ManifoldConf>::update(this->getControlVertex(), dpV);
			};

			/// \brief Revert the last state update.
			virtual void revertUpdateImplementation() { this->getControlVertex() = _p_v ; }

			inline aslam::backend::DesignVariable & getDesignVariable(){ return *this; }
			inline const aslam::backend::DesignVariable & getDesignVariable() const { return *this; }
		};
	}
}

namespace aslam {
namespace splines {

template <typename TBSplineConf, typename TDerivedConf = TBSplineConf>
class OPTBSpline : public internal::DesignVariableSegmentBSplineConf<TDerivedConf>::BSpline {
protected:
	typedef internal::DesignVariableSegmentBSplineConf<TDerivedConf> configuration_t;
public:
	typedef typename configuration_t::BSpline parent_t;
	typedef OPTBSpline<TDerivedConf> spline_t;
	typedef typename parent_t::TimePolicy TimePolicy;
	typedef typename parent_t::time_t time_t;
	typedef typename parent_t::point_t point_t;
	typedef typename parent_t::SegmentIterator SegmentIterator;
	typedef typename parent_t::SegmentConstIterator SegmentConstIterator;

	typedef TDerivedConf CONF;

	typedef spline_t TYPE;

	enum {
		PointSize = TBSplineConf::PointSize::VALUE,
		Dimension = TBSplineConf::Dimension::VALUE
	};

	typedef aslam::backend::DesignVariable dv_t;
	typedef aslam::backend::VectorExpression<PointSize> expression_t;


	OPTBSpline(const CONF & config) : parent_t(configuration_t(config)){}

	void init();

	size_t numDesignVariables();
	dv_t * designVariable(size_t i);

	std::vector<dv_t *> getDesignVariables(time_t time);

	// add one Segment at the end of the PoseSpline
	void addSegment(time_t t, const point_t & p);
	void addSegment2(time_t t, const point_t & p, double lambda);
	void removeSegment();

	template <int IMaxDerivativeOrder>
	class ExpressionFactory {
	protected:
		typedef typename spline_t::template Evaluator<IMaxDerivativeOrder> eval_t;
		typedef boost::shared_ptr<const eval_t> eval_ptr_t;
		eval_ptr_t _evalPtr;

	public:
		ExpressionFactory(spline_t & spline, const time_t & t);

		const eval_t & getEvaluator() const { return *_evalPtr; };

		/// \brief get an expression
		expression_t toExpression(int derivativeOrder);
	};

	template <int IMaxDerivativeOrder> inline ExpressionFactory<IMaxDerivativeOrder> getExpressionFactoryAt(const time_t & t){ return ExpressionFactory<IMaxDerivativeOrder>(*this, t); }
protected:
	/// \brief the vector of design variables.
	std::vector< dv_t * > _designVariables;
};

} // namespace splines
} // namespace aslam

#include "implementation/OPTBSplineImpl.hpp"

#endif /* OPTBSPLINE_HPP_ */
